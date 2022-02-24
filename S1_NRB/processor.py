import os
import re
import time
import tempfile
from getpass import getpass
from datetime import datetime
from lxml import etree
import numpy as np
from osgeo import gdal
from spatialist import Raster, Vector, vectorize, boundary, bbox, intersect, rasterize
from spatialist.ancillary import finder
from spatialist.auxil import gdalwarp, gdalbuildvrt
from pyroSAR import identify_many, Archive
from pyroSAR.snap.util import geocode, noise_power
from pyroSAR.ancillary import groupbyTime, seconds, find_datasets
from pyroSAR.auxdata import dem_autoload, dem_create
from S1_NRB.config import get_config, geocode_conf, gdal_conf
import S1_NRB.ancillary as ancil
import S1_NRB.tile_extraction as tile_ex
from S1_NRB.metadata import extract, xmlparser, stacparser
gdal.UseExceptions()


def nrb_processing(config, scenes, datadir, outdir, tile, extent, epsg, wbm=None, multithread=True,
                   compress='LERC_ZSTD', overviews=None, recursive=False):
    """
    Finalizes the generation of Sentinel-1 NRB products after the main processing steps via `pyroSAR.snap.util.geocode`
    have been executed. This includes the following:
    - Converting all GeoTIFF files to Cloud Optimized GeoTIFF (COG) format
    - Generating annotation datasets in Virtual Raster Tile (VRT) format
    - Applying the NRB product directory structure & naming convention
    - Generating metadata in XML and STAC JSON formats
    
    Parameters
    ----------
    config: dict
        Dictionary of the parsed config parameters for the current process.
    scenes: list[str]
        List of scenes to process. Either an individual scene or multiple, matching scenes (consecutive acquisitions).
    datadir: str
        The directory containing the datasets processed from the source scenes using pyroSAR.
    outdir: str
        The directory to write the final files to.
    tile: str
        ID of an MGRS tile.
    extent: dict
        Spatial extent of the MGRS tile, derived from a `spatialist.vector.Vector` object.
    epsg: int
        The CRS used for the NRB product; provided as an EPSG code.
    wbm: str, optional
        Path to a water body mask file with the dimensions of an MGRS tile.
    multithread: bool, optional
        Should `gdalwarp` use multithreading? Default is True. The number of threads used, can be adjusted in the
        config.ini file with the parameter `gdal_threads`.
    compress: str, optional
        Compression algorithm to use. See https://gdal.org/drivers/raster/gtiff.html#creation-options for options.
        Defaults to 'LERC_DEFLATE'.
    overviews: list[int], optional
        Internal overview levels to be created for each GeoTIFF file. Defaults to [2, 4, 9, 18, 36]
    recursive: bool, optional
        Find datasets recursively in the directory specified with the parameter `datadir`? Default is False.
    
    Returns
    -------
    None
    """
    if overviews is None:
        overviews = [2, 4, 9, 18, 36]
    
    ids = identify_many(scenes)
    datasets = [find_datasets(directory=datadir,
                              recursive=recursive,
                              start=id.start, stop=id.start)
                for id in ids]
    
    if len(datasets) == 0:
        raise RuntimeError("No pyroSAR datasets were found in the directory '{}'".format(datadir))
    
    pattern = '[VH]{2}_gamma0-rtc'
    i = 0
    snap_dm_tile_overlap = []
    while i < len(datasets):
        pols = [x for x in datasets[i] if re.search(pattern, os.path.basename(x))]
        snap_dm_ras = re.sub(pattern, 'datamask', pols[0])
        snap_dm_vec = snap_dm_ras.replace('.tif', '.gpkg')
        
        if not all([os.path.isfile(x) for x in [snap_dm_ras, snap_dm_vec]]):
            with Raster(pols[0]) as ras:
                arr = ras.array()
                mask = ~np.isnan(arr)
                with vectorize(target=mask, reference=ras) as vec:
                    with boundary(vec, expression="value=1") as bounds:
                        if not os.path.isfile(snap_dm_ras):
                            print('creating raster mask', i)
                            rasterize(vectorobject=bounds, reference=ras, outname=snap_dm_ras)
                        if not os.path.isfile(snap_dm_vec):
                            print('creating vector mask', i)
                            bounds.write(outfile=snap_dm_vec)
        with Vector(snap_dm_vec) as bounds:
            with bbox(extent, epsg) as tile_geom:
                inter = intersect(bounds, tile_geom)
                if inter is None:
                    print('removing dataset', i)
                    del ids[i]
                    del datasets[i]
                else:
                    # Add snap_dm_ras to list if it overlaps with the current tile
                    snap_dm_tile_overlap.append(snap_dm_ras)
                    i += 1
                    inter.close()
    
    if len(ids) == 0:
        raise RuntimeError('None of the scenes overlap with the current tile {tile_id}: '
                           '\n{scenes}'.format(tile_id=tile, scenes=scenes))
    
    src_scenes = [i.scene for i in ids]
    product_start, product_stop = ancil.calc_product_start_stop(src_scenes=src_scenes, extent=extent, epsg=epsg)
    
    meta = {'mission': ids[0].sensor,
            'mode': ids[0].meta['acquisition_mode'],
            'start': product_start,
            'orbitnumber': ids[0].meta['orbitNumbers_abs']['start'],
            'datatake': hex(ids[0].meta['frameNumber']).replace('x', '').upper(),
            'stop': product_stop,
            'tile': tile,
            'id': 'ABCD'}
    
    skeleton = '{mission}_{mode}_NRB__1SDV_{start}_{stop}_{orbitnumber:06}_{datatake}_{tile}_{id}'
    
    nrbdir = os.path.join(outdir, skeleton.format(**meta))
    os.makedirs(nrbdir, exist_ok=True)
    
    metaL = meta.copy()
    for key, val in metaL.items():
        if not isinstance(val, int):
            metaL[key] = val.lower()
    
    skeleton = '{mission}-{mode}-nrb-{start}-{stop}-{orbitnumber:06}-{datatake}-{tile}-{suffix}.tif'
    
    # 'z_error': Maximum error threshold on values for LERC* compression.
    # Will be ignored if a compression algorithm is used that isn't related to LERC.
    item_map = {'VV_gamma0': {'suffix': 'vv-g-lin',
                              'z_error': 1e-4},
                'VH_gamma0': {'suffix': 'vh-g-lin',
                              'z_error': 1e-4},
                'HH_gamma0': {'suffix': 'hh-g-lin',
                              'z_error': 1e-4},
                'HV_gamma0': {'suffix': 'hv-g-lin',
                              'z_error': 1e-4},
                'incidenceAngleFromEllipsoid': {'suffix': 'ei',
                                                'z_error': 1e-3},
                'layoverShadowMask': {'suffix': 'dm',
                                      'z_error': 0},
                'localIncidenceAngle': {'suffix': 'li',
                                        'z_error': 1e-2},
                'scatteringArea': {'suffix': 'lc',
                                   'z_error': 0.1},
                'gammaSigmaRatio': {'suffix': 'gs',
                                    'z_error': 1e-4},
                'acquisitionImage': {'suffix': 'id',
                                     'z_error': 0},
                'VV_NESZ': {'suffix': 'np-vv',
                            'z_error': 2e-5},
                'VH_NESZ': {'suffix': 'np-vh',
                            'z_error': 2e-5},
                'HH_NESZ': {'suffix': 'np-hh',
                            'z_error': 2e-5},
                'HV_NESZ': {'suffix': 'np-hv',
                            'z_error': 2e-5}}
    
    driver = 'COG'
    ovr_resampling = 'AVERAGE'
    write_options_base = ['BLOCKSIZE=512', 'OVERVIEW_RESAMPLING={}'.format(ovr_resampling)]
    write_options = dict()
    for key in item_map:
        write_options[key] = write_options_base.copy()
        if compress is not None:
            entry = 'COMPRESS={}'.format(compress)
            write_options[key].append(entry)
            if compress.startswith('LERC'):
                entry = 'MAX_Z_ERROR={:f}'.format(item_map[key]['z_error'])
                write_options[key].append(entry)
    
    ####################################################################################################################
    # format existing datasets found by `pyroSAR.ancillary.find_datasets`
    if len(datasets) > 1:
        files = list(zip(*datasets))
    else:
        files = datasets[0]
    
    pattern = '|'.join(item_map.keys())
    for i, item in enumerate(files):
        if isinstance(item, str):
            match = re.search(pattern, item)
            if match is not None:
                key = match.group()
            else:
                continue
        else:
            match = [re.search(pattern, x) for x in item]
            keys = [x if x is None else x.group() for x in match]
            if len(list(set(keys))) != 1:
                raise RuntimeError('file mismatch:\n{}'.format('\n'.join(item)))
            if None in keys:
                continue
            key = keys[0]
        
        if key == 'layoverShadowMask':
            # The data mask will be created later on in the processing workflow.
            continue
        
        metaL['suffix'] = item_map[key]['suffix']
        outname_base = skeleton.format(**metaL)
        if re.search('_gamma0', key):
            subdir = 'measurement'
        else:
            subdir = 'annotation'
        outname = os.path.join(nrbdir, subdir, outname_base)
        
        if not os.path.isfile(outname):
            os.makedirs(os.path.dirname(outname), exist_ok=True)
            print(outname)
            bounds = [extent['xmin'], extent['ymin'], extent['xmax'], extent['ymax']]
            
            if isinstance(item, tuple):
                with Raster(list(item), list_separate=False) as ras:
                    source = ras.filename
            elif isinstance(item, str):
                source = tempfile.NamedTemporaryFile(suffix='.vrt').name
                gdalbuildvrt(item, source)
            else:
                raise TypeError('type {} is not supported: {}'.format(type(item), item))
            
            # modify temporary VRT to make sure overview levels and resampling are properly applied
            tree = etree.parse(source)
            root = tree.getroot()
            ovr = etree.SubElement(root, 'OverviewList', attrib={'resampling': ovr_resampling.lower()})
            ov = str(overviews)
            for x in ['[', ']', ',']:
                ov = ov.replace(x, '')
            ovr.text = ov
            etree.indent(root)
            tree.write(source, pretty_print=True, xml_declaration=False, encoding='utf-8')
            
            snap_nodata = 0
            gdalwarp(source, outname,
                     options={'format': driver, 'outputBounds': bounds, 'srcNodata': snap_nodata, 'dstNodata': 'nan',
                              'multithread': multithread, 'creationOptions': write_options[key]})
    
    proc_time = datetime.now()
    t = proc_time.isoformat().encode()
    product_id = ancil.generate_unique_id(encoded_str=t)
    nrbdir_new = nrbdir.replace('ABCD', product_id)
    os.rename(nrbdir, nrbdir_new)
    nrbdir = nrbdir_new
    
    if type(files[0]) == tuple:
        files = [item for tup in files for item in tup]
    gs_path = finder(nrbdir, [r'gs\.tif$'], regex=True)[0]
    measure_paths = finder(nrbdir, ['[hv]{2}-g-lin.tif$'], regex=True)
    
    ####################################################################################################################
    # log-scaled gamma nought & color composite
    log_vrts = []
    for item in measure_paths:
        log = item.replace('lin.tif', 'log.vrt')
        log_vrts.append(log)
        if not os.path.isfile(log):
            print(log)
            ancil.vrt_pixfun(src=item, dst=log, fun='log10', scale=10,
                             options={'VRTNodata': 'NaN'}, overviews=overviews, overview_resampling=ovr_resampling)
    
    cc_path = re.sub('[hv]{2}', 'cc', measure_paths[0]).replace('.tif', '.vrt')
    # cc_path = re.sub('[hv]{2}', 'cc', log_vrts[0])
    ancil.create_rgb_vrt(outname=cc_path, infiles=measure_paths, overviews=overviews,
                         overview_resampling=ovr_resampling)
    
    ####################################################################################################################
    # Data mask
    if not config['dem_type'] == 'GETASSE30' and not os.path.isfile(wbm):
        raise FileNotFoundError('External water body mask could not be found: {}'.format(wbm))
    
    dm_path = gs_path.replace('-gs.tif', '-dm.tif')
    ancil.create_data_mask(outname=dm_path, valid_mask_list=snap_dm_tile_overlap, src_files=files,
                           extent=extent, epsg=epsg, driver=driver, creation_opt=write_options['layoverShadowMask'],
                           overviews=overviews, overview_resampling=ovr_resampling, wbm=wbm)
    
    ####################################################################################################################
    # Acquisition ID image
    ancil.create_acq_id_image(ref_tif=gs_path, valid_mask_list=snap_dm_tile_overlap, src_scenes=src_scenes,
                              extent=extent, epsg=epsg, driver=driver, creation_opt=write_options['acquisitionImage'],
                              overviews=overviews)
    
    ####################################################################################################################
    # sigma nought RTC
    for item in measure_paths:
        sigma0_rtc_lin = item.replace('g-lin.tif', 's-lin.vrt')
        sigma0_rtc_log = item.replace('g-lin.tif', 's-log.vrt')
        
        if not os.path.isfile(sigma0_rtc_lin):
            print(sigma0_rtc_lin)
            ancil.vrt_pixfun(src=[item, gs_path], dst=sigma0_rtc_lin, fun='mul',
                             options={'VRTNodata': 'NaN'}, overviews=overviews, overview_resampling=ovr_resampling)
            ancil.vrt_relpath(sigma0_rtc_lin)
        
        if not os.path.isfile(sigma0_rtc_log):
            print(sigma0_rtc_log)
            ancil.vrt_pixfun(src=sigma0_rtc_lin, dst=sigma0_rtc_log, fun='log10', scale=10,
                             options={'VRTNodata': 'NaN'}, overviews=overviews, overview_resampling=ovr_resampling)
    
    ####################################################################################################################
    # metadata
    nrb_tifs = finder(nrbdir, ['-[a-z]{2,3}.tif'], regex=True, recursive=True)
    meta = extract.meta_dict(config=config, target=nrbdir, src_scenes=src_scenes, src_files=files, proc_time=proc_time)
    xmlparser.main(meta=meta, target=nrbdir, tifs=nrb_tifs)
    stacparser.main(meta=meta, target=nrbdir, tifs=nrb_tifs)


def prepare_dem(id_list, config, threads, epsg, tr, buffer=None):
    """
    Downloads and prepares the chosen DEM for the NRB processing workflow.
    
    Parameters
    ----------
    id_list: list[pyroSAR.drivers.ID]
        A list of pyroSAR metadata handlers for the current selection of scenes.
    config: dict
        Dictionary of the parsed config parameters for the current process.
    threads: int
        The number of threads to pass to `pyroSAR.auxdata.dem_create`.
    epsg: int
        The CRS used for the NRB product; provided as an EPSG code.
    tr: int
        The target resolution to pass to `pyroSAR.auxdata.dem_create`.
    buffer: float, optional
        A buffer in degrees to pass to `pyroSAR.auxdata.dem_autoload`.
    
    Returns
    -------
    dem_names: list[list[str]]
        List of lists containing paths to DEM tiles for each scene to be processed.
    """
    if config['dem_type'] == 'GETASSE30':
        geoid = 'WGS84'
    else:
        geoid = 'EGM2008'
    
    wbm_dems = ['Copernicus 10m EEA DEM',
                'Copernicus 30m Global DEM II']
    wbm_dir = os.path.join(config['wbm_dir'], config['dem_type'])
    dem_dir = os.path.join(config['dem_dir'], config['dem_type'])
    username = None
    password = None
    
    boxes = [i.bbox() for i in id_list]
    max_ext = ancil.get_max_ext(boxes=boxes, buffer=buffer)
    ext_id = ancil.generate_unique_id(encoded_str=str(max_ext).encode())
    
    fname_wbm_tmp = os.path.join(wbm_dir, 'mosaic_{}.vrt'.format(ext_id))
    fname_dem_tmp = os.path.join(dem_dir, 'mosaic_{}.vrt'.format(ext_id))
    if config['dem_type'] in wbm_dems:
        if not os.path.isfile(fname_wbm_tmp) or not os.path.isfile(fname_dem_tmp):
            username = input('Please enter your DEM access username:')
            password = getpass('Please enter your DEM access password:')
        os.makedirs(wbm_dir, exist_ok=True)
        if not os.path.isfile(fname_wbm_tmp):
            dem_autoload(boxes, demType=config['dem_type'],
                         vrt=fname_wbm_tmp, buffer=buffer, product='wbm',
                         username=username, password=password,
                         nodata=1, hide_nodata=True)
    os.makedirs(dem_dir, exist_ok=True)
    if not os.path.isfile(fname_dem_tmp):
        dem_autoload(boxes, demType=config['dem_type'],
                     vrt=fname_dem_tmp, buffer=buffer, product='dem',
                     username=username, password=password)
    
    dem_names = []
    for scene in id_list:
        dem_names_scene = []
        with scene.bbox() as box:
            tiles = tile_ex.tiles_from_aoi(vectorobject=box, kml=config['kml_file'], epsg=epsg)
            for tilename in tiles:
                dem_tile = os.path.join(dem_dir, '{}_DEM.tif'.format(tilename))
                dem_names_scene.append(dem_tile)
                if not os.path.isfile(dem_tile):
                    with tile_ex.extract_tile(config['kml_file'], tilename) as tile:
                        ext = tile.extent
                        bounds = [ext['xmin'], ext['ymin'], ext['xmax'], ext['ymax']]
                        dem_create(src=fname_dem_tmp, dst=dem_tile, t_srs=epsg, tr=(tr, tr),
                                   geoid_convert=True, geoid=geoid, pbar=True,
                                   outputBounds=bounds, threads=threads)
            if os.path.isfile(fname_wbm_tmp):
                for tilename in tiles:
                    wbm_tile = os.path.join(wbm_dir, '{}_WBM.tif'.format(tilename))
                    if not os.path.isfile(wbm_tile):
                        with tile_ex.extract_tile(config['kml_file'], tilename) as tile:
                            ext = tile.extent
                            bounds = [ext['xmin'], ext['ymin'], ext['xmax'], ext['ymax']]
                            dem_create(src=fname_wbm_tmp, dst=wbm_tile, t_srs=epsg, tr=(tr, tr),
                                       resampling_method='mode', pbar=True,
                                       outputBounds=bounds, threads=threads)
        dem_names.append(dem_names_scene)
    boxes = None  # make sure all bounding box Vector objects are deleted
    
    return dem_names


def main(config_file, section_name):
    config = get_config(config_file=config_file, section_name=section_name)
    log = ancil.set_logging(config=config)
    geocode_prms = geocode_conf(config=config)
    gdal_prms = gdal_conf(config=config)
    
    snap_flag = True
    nrb_flag = True
    if config['mode'] == 'snap':
        nrb_flag = False
    elif config['mode'] == 'nrb':
        snap_flag = False
    
    ####################################################################################################################
    # archive / scene selection
    scenes = finder(config['scene_dir'], [r'^S1[AB].*\.zip'], regex=True, recursive=True)
    if not os.path.isfile(config['db_file']):
        config['db_file'] = os.path.join(config['work_dir'], config['db_file'])
    
    with Archive(dbfile=config['db_file']) as archive:
        archive.insert(scenes)
        selection = archive.select(product='SLC',
                                   acquisition_mode=config['acq_mode'],
                                   mindate=config['mindate'], maxdate=config['maxdate'])
    ids = identify_many(selection)
    
    if len(selection) == 0:
        raise RuntimeError("No scenes could be found for acquisition mode '{acq_mode}', mindate '{mindate}' "
                           "and maxdate '{maxdate}' in directory '{scene_dir}'.".format(acq_mode=config['acq_mode'],
                                                                                        mindate=config['mindate'],
                                                                                        maxdate=config['maxdate'],
                                                                                        scene_dir=config['scene_dir']))
    
    ####################################################################################################################
    # geometry and DEM handling
    geo_dict, align_dict = tile_ex.main(config=config, tr=geocode_prms['tr'])
    aoi_tiles = list(geo_dict.keys())
    
    epsg_set = set([geo_dict[tile]['epsg'] for tile in list(geo_dict.keys())])
    if len(epsg_set) != 1:
        raise RuntimeError('The AOI covers multiple UTM zones: {}\n '
                           'This is currently not supported. Please refine your AOI.'.format(list(epsg_set)))
    epsg = epsg_set.pop()
    
    dem_names = prepare_dem(id_list=ids, config=config, threads=gdal_prms['threads'], epsg=epsg,
                            tr=geocode_prms['tr'], buffer=1.5)
    
    ####################################################################################################################
    # geocode & noise power - SNAP processing
    np_dict = {'sigma0': 'NESZ', 'beta0': 'NEBZ', 'gamma0': 'NEGZ'}
    np_refarea = 'sigma0'
    
    if snap_flag:
        for i, scene in enumerate(ids):
            dem_buffer = 200  # meters
            fname_dem = os.path.join(config['tmp_dir'], scene.outname_base() + '_DEM.tif')
            if not os.path.isfile(fname_dem):
                with scene.geometry() as footprint:
                    footprint.reproject(epsg)
                    extent = footprint.extent
                    extent['xmin'] -= dem_buffer
                    extent['ymin'] -= dem_buffer
                    extent['xmax'] += dem_buffer
                    extent['ymax'] += dem_buffer
                    with bbox(extent, epsg) as dem_box:
                        with Raster(dem_names[i], list_separate=False)[dem_box] as dem_mosaic:
                            dem_mosaic.write(fname_dem, format='GTiff')
            
            list_processed = finder(config['out_dir'], [scene.start], regex=True, recursive=False)
            exclude = list(np_dict.values())
            print('###### [GEOCODE] Scene {s}/{s_total}: {scene}'.format(s=i+1, s_total=len(ids),
                                                                         scene=scene.scene))
            if len([item for item in list_processed if not any(ex in item for ex in exclude)]) < 3:
                start_time = time.time()
                try:
                    geocode(infile=scene, outdir=config['out_dir'], t_srs=epsg, tmpdir=config['tmp_dir'],
                            standardGridOriginX=align_dict['xmax'], standardGridOriginY=align_dict['ymin'],
                            externalDEMFile=fname_dem, **geocode_prms)
                    
                    t = round((time.time() - start_time), 2)
                    log.info('[GEOCODE] -- {scene} -- {time}'.format(scene=scene.scene, time=t))
                    if t <= 500:
                        log.warning('[GEOCODE] -- {scene} -- Processing might have terminated prematurely. Check'
                                    ' terminal for uncaught SNAP errors!'.format(scene=scene.scene))
                except Exception as e:
                    log.error('[GEOCODE] -- {scene} -- {error}'.format(scene=scene.scene, error=e))
                    continue
            else:
                msg = 'Already processed - Skip!'
                print('### ' + msg)
                log.info('[GEOCODE] -- {scene} -- {msg}'.format(scene=scene.scene, msg=msg))
            
            print('###### [NOISE_P] Scene {s}/{s_total}: {scene}'.format(s=i+1, s_total=len(ids),
                                                                         scene=scene.scene))
            if len([item for item in list_processed if np_dict[np_refarea] in item]) == 0:
                start_time = time.time()
                try:
                    noise_power(infile=scene.scene, outdir=config['out_dir'], polarizations=scene.polarizations,
                                spacing=geocode_prms['tr'], t_srs=epsg, refarea=np_refarea, tmpdir=config['tmp_dir'],
                                externalDEMFile=fname_dem, externalDEMApplyEGM=geocode_prms['externalDEMApplyEGM'],
                                alignToStandardGrid=geocode_prms['alignToStandardGrid'],
                                standardGridOriginX=align_dict['xmax'], standardGridOriginY=align_dict['ymin'],
                                clean_edges=geocode_prms['clean_edges'],
                                clean_edges_npixels=geocode_prms['clean_edges_npixels'],
                                rlks=geocode_prms['rlks'], azlks=geocode_prms['azlks'])
                    t = round((time.time() - start_time), 2)
                    log.info('[NOISE_P] -- {scene} -- {time}'.format(scene=scene.scene, time=t))
                except Exception as e:
                    log.error('[NOISE_P] -- {scene} -- {error}'.format(scene=scene.scene, error=e))
                    continue
            else:
                msg = 'Already processed - Skip!'
                print('### ' + msg)
                log.info('[NOISE_P] -- {scene} -- {msg}'.format(scene=scene.scene, msg=msg))
    
    ####################################################################################################################
    # NRB - final product generation
    if nrb_flag:
        selection_grouped = groupbyTime(images=selection, function=seconds, time=60)
        for t, tile in enumerate(aoi_tiles):
            outdir = os.path.join(config['out_dir'], tile)
            os.makedirs(outdir, exist_ok=True)
            wbm = os.path.join(config['wbm_dir'], config['dem_type'], '{}_WBM.tif'.format(tile))
            if not os.path.isfile(wbm):
                wbm = None
            
            for s, scenes in enumerate(selection_grouped):
                if isinstance(scenes, str):
                    scenes = [scenes]
                print('###### [NRB] Tile {t}/{t_total}: {tile} | '
                      'Scenes {s}/{s_total}: {scenes} '.format(tile=tile, t=t+1, t_total=len(aoi_tiles),
                                                               scenes=[os.path.basename(s) for s in scenes],
                                                               s=s+1, s_total=len(selection_grouped)))
                start_time = time.time()
                try:
                    nrb_processing(config=config, scenes=scenes, datadir=os.path.dirname(outdir), outdir=outdir,
                                   tile=tile, extent=geo_dict[tile]['ext'], epsg=epsg, wbm=wbm,
                                   multithread=gdal_prms['multithread'])
                    log.info('[    NRB] -- {scenes} -- {time}'.format(scenes=scenes,
                                                                      time=round((time.time() - start_time), 2)))
                except Exception as e:
                    log.exception('[    NRB] -- {scenes} -- {error}'.format(scenes=scenes, error=e))
                    continue
        
        gdal.SetConfigOption('GDAL_NUM_THREADS', gdal_prms['threads_before'])
