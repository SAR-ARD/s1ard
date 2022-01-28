import os
import re
import time
import numpy as np
from osgeo import gdal
from spatialist import Raster, Vector, vectorize, boundary, bbox, intersect, rasterize
from spatialist.ancillary import finder
from spatialist.auxil import gdalwarp
from pyroSAR import identify_many, Archive
from pyroSAR.snap.util import geocode, noise_power
from pyroSAR.ancillary import groupbyTime, seconds, find_datasets

from S1_NRB.config import get_config, geocode_params
import S1_NRB.ancillary as ancil
import S1_NRB.tile_extraction as tile_ex
from S1_NRB.metadata import extract, xmlparser, stacparser

gdal.UseExceptions()


def nrb_processing(scenes, datadir, outdir, tile, extent, epsg, dem_name, compress='LERC_ZSTD',
                   overviews=None, recursive=False, external_wbm=None):
    """
    Finalizes the generation of Sentinel-1 NRB products after the main processing steps via `pyroSAR.snap.util.geocode`
    have been executed. This includes the following:
    - Converting all GeoTIFF files to Cloud Optimized GeoTIFF (COG) format
    - Generating annotation datasets in Virtual Raster Tile (VRT) format
    - Applying the NRB product directory structure & naming convention
    - Generating metadata in XML and STAC JSON formats
    
    Parameters
    ----------
    scenes: list[str]
        List of scenes to process. Either an individual scene or multiple, matching scenes (consecutive acquisitions).
    datadir: str
        the directory containing the datasets processed from the scenes using pyroSAR.
    outdir: str
        The directory to write the final files to.
    tile: str
        ID of an MGRS tile.
    extent: dict
        Spatial extent of the MGRS tile, derived from a `spatialist.vector.Vector` object.
    epsg: int
        The CRS used for the NRB product; provided as an EPSG code.
    dem_name: str
        Name of the DEM used for processing. Must match with supported options listed by `pyroSAR.snap.util.geocode`
    compress: str
        Compression algorithm to use. See https://gdal.org/drivers/raster/gtiff.html#creation-options for options.
        Defaults to 'LERC_DEFLATE'.
    overviews: list[int], optional
        Internal overview levels to be created for each GeoTIFF file. Defaults to [2, 4, 8, 16, 32]
    recursive: bool, optional
        Find datasets recursively in the directory specified with the parameter 'outdir'? Default is False.
    external_wbm: str, optional
        Path to the COP-DEM water body mask.
    
    Returns
    -------
    None
    """
    
    if overviews is None:
        overviews = [2, 4, 8, 16, 32]
    
    wbm = False
    if external_wbm is not None:
        wbm = True
        if not os.path.isfile(external_wbm):
            raise FileNotFoundError('External water body mask could not be found: {}'.format(external_wbm))
    
    ids = identify_many(scenes)
    
    datasets = [find_datasets(directory=datadir,
                              recursive=recursive,
                              start=id.start, stop=id.start)
                for id in ids]
    
    if len(datasets) == 0:
        raise RuntimeError("No pyroSAR datasets were found in the directory '{}'".format(datadir))
    
    pattern = '[VH]{2}_gamma0-rtc'
    
    i = 0
    datamask_ras_list = []
    while i < len(datasets):
        pols = [x for x in datasets[i] if re.search(pattern, os.path.basename(x))]
        datamask_ras = re.sub(pattern, 'datamask', pols[0])
        datamask_vec = datamask_ras.replace('.tif', '.gpkg')
        
        if not all([os.path.isfile(x) for x in [datamask_ras, datamask_vec]]):
            with Raster(pols[0]) as ras:
                arr = ras.array()
                mask = ~np.isnan(arr)
                with vectorize(target=mask, reference=ras) as vec:
                    with boundary(vec, expression="value=1") as bounds:
                        if not os.path.isfile(datamask_ras):
                            print('creating raster mask', i)
                            rasterize(vectorobject=bounds, reference=ras, outname=datamask_ras)
                        if not os.path.isfile(datamask_vec):
                            print('creating vector mask', i)
                            bounds.write(outfile=datamask_vec)
        
        with Vector(datamask_vec) as bounds:
            with bbox(extent, epsg) as tile_geom:
                inter = intersect(bounds, tile_geom)
                if inter is None:
                    print('removing dataset', i)
                    del ids[i]
                    del datasets[i]
                else:
                    datamask_ras_list.append(datamask_ras)
                    i += 1
                    inter.close()
    
    if len(ids) == 0:
        raise RuntimeError('none of the scenes overlaps with the tile')
    
    starts = [id.start for id in ids]
    stops = [id.stop for id in ids]
    
    meta = {'mission': ids[0].sensor,
            'mode': ids[0].meta['acquisition_mode'],
            'start': min(starts),
            'orbitnumber': ids[0].meta['orbitNumbers_abs']['start'],
            'datatake': hex(ids[0].meta['frameNumber']).replace('x', '').upper(),
            'stop': max(stops),
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
                'VV_NEGZ': {'suffix': 'np-vv',
                            'z_error': 2e-5},
                'VH_NEGZ': {'suffix': 'np-vh',
                            'z_error': 2e-5},
                'HH_NEGZ': {'suffix': 'np-hh',
                            'z_error': 2e-5},
                'HV_NEGZ': {'suffix': 'np-hv',
                            'z_error': 2e-5}}
    
    driver = 'COG'
    write_options_base = ['BLOCKSIZE=512', 'OVERVIEW_RESAMPLING=AVERAGE']
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
        
        val = item_map[key]['suffix']
        metaL['suffix'] = val
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
            
            if isinstance(item, str):
                source = item
            else:
                # mosaic files in temporary VRT
                with Raster(list(item), list_separate=False) as ras:
                    source = ras.filename
            
            if outname.endswith('-dm.tif'):
                dstnodata = 255
            else:
                dstnodata = 'nan'
            
            gdalwarp(source, outname,
                     options={'format': driver, 'outputBounds': bounds, 'srcNodata': 0, 'dstNodata': dstnodata,
                              'creationOptions': write_options[key]})
    
    product_id, proc_time = ancil.generate_product_id()
    nrbdir_new = nrbdir.replace('ABCD', product_id)
    os.rename(nrbdir, nrbdir_new)
    nrbdir = nrbdir_new
    
    if type(files[0]) == tuple:
        src_files = [item for tup in files for item in tup]
    else:
        src_files = files
    
    ####################################################################################################################
    # log-scaled gamma nought
    
    measure = finder(nrbdir, ['[hv]{2}-g-lin.tif$'], regex=True)
    for item in measure:
        log = item.replace('lin.tif', 'log.vrt')
        if not os.path.isfile(log):
            print(log)
            ancil.vrt_pixfun(src=item,
                             dst=log,
                             fun='log10',
                             scale=10,
                             options={'VRTNodata': 'NaN'})
    
    ####################################################################################################################
    # modify data mask
    
    dm_path = finder(nrbdir, [r'dm\.tif$'], regex=True)[0]
    ancil.modify_data_mask(dm_path=dm_path, mask_list=datamask_ras_list, src_files=src_files, extent=extent,
                           epsg=epsg, driver=driver, creation_opt=write_options['layoverShadowMask'],
                           overviews=overviews, multilayer=True, wbm=wbm, wbm_path=external_wbm)
    
    ####################################################################################################################
    # Acquisition ID image
    ancil.create_acq_id_image(outdir=os.path.join(nrbdir, 'annotation'), ref_tif=dm_path, src_scenes=scenes,
                              src_files=src_files, driver=driver, creation_opt=write_options['acquisitionImage'],
                              overviews=overviews)
    
    ####################################################################################################################
    # sigma nought RTC
    
    gs_name = finder(nrbdir, [r'gs\.tif$'], regex=True)[0]
    for item in measure:
        sigma0_rtc_lin = item.replace('g-lin.tif', 's-lin.vrt')
        sigma0_rtc_log = item.replace('g-lin.tif', 's-log.vrt')
        
        if not os.path.isfile(sigma0_rtc_lin):
            print(sigma0_rtc_lin)
            ancil.vrt_pixfun(src=[item, gs_name],
                             dst=sigma0_rtc_lin,
                             fun='mul',
                             options={'VRTNodata': 'NaN'})
            ancil.vrt_relpath(sigma0_rtc_lin)
        
        if not os.path.isfile(sigma0_rtc_log):
            print(sigma0_rtc_log)
            ancil.vrt_pixfun(src=sigma0_rtc_lin,
                             dst=sigma0_rtc_log,
                             fun='log10',
                             scale=10,
                             options={'VRTNodata': 'NaN'})
    
    ####################################################################################################################
    # metadata
    tifs = finder(nrbdir, ['-[a-z]{2,3}.tif'], regex=True)
    meta = extract.meta_dict(target=nrbdir, src_scenes=scenes, src_files=src_files,
                             dem_name=dem_name, proc_time=proc_time)
    xmlparser.main(meta=meta, target=nrbdir, tifs=tifs)
    stacparser.main(meta=meta, target=nrbdir, tifs=tifs)


def main(config_file, section_name):
    config = get_config(config_file=config_file, section_name=section_name)
    log = ancil.set_logging(config=config)
    geocode_prms = geocode_params(config=config)
    
    geocode_flag = True
    nrb_flag = True
    if config['mode'] == 'snap':
        nrb_flag = False
    elif config['mode'] == 'nrb':
        geocode_flag = False
    
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
    
    # avoid reprocessing of already (successfully) SNAP processed scenes when modes 'all' or 'snap' are used
    if config['mode'] != 'nrb':
        selection = ancil.filter_selection(selection=selection, processdir=config['out_dir'])
    
    if len(selection) == 0:
        raise RuntimeError("No scenes could be found for acquisition mode '{acq_mode}', mindate '{mindate}' "
                           "and maxdate '{maxdate}' in directory '{scene_dir}'.".format(acq_mode=config['acq_mode'],
                                                                                        mindate=config['mindate'],
                                                                                        maxdate=config['maxdate'],
                                                                                        scene_dir=config['scene_dir']))
    ####################################################################################################################
    # geometry handling
    
    geo_dict, align_dict = tile_ex.main(config=config, tr=geocode_prms['tr'])
    tiles = list(geo_dict.keys())
    
    epsg_set = set([geo_dict[tile]['epsg'] for tile in list(geo_dict.keys())])
    if len(epsg_set) != 1:
        raise RuntimeError('The AOI covers multiple UTM zones: {}\n '
                           'This is currently not supported. Please refine your AOI.'.format(list(epsg_set)))
    epsg = epsg_set.pop()
    
    ####################################################################################################################
    # geocode & noise power - SNAP processing
    
    if geocode_flag:
        # Process scenes individually
        for scene in selection:
            print('###### SNAP GEOCODE: {scene}'.format(scene=scene))
            start_time = time.time()
            try:
                geocode(infile=scene, outdir=config['out_dir'], t_srs=epsg, tmpdir=config['tmp_dir'],
                        standardGridOriginX=align_dict['xmax'], standardGridOriginY=align_dict['ymin'],
                        externalDEMFile=config.get('ext_dem_file'), **geocode_prms)
                
                t = round((time.time() - start_time), 2)
                log.info('[GEOCODE] -- {scene} -- {time}'.format(scene=scene, time=t))
                if t <= 500:
                    log.warning('[GEOCODE] -- {scene} -- Processing might have terminated prematurely. Check'
                                ' terminal for uncaught SNAP errors!'.format(scene=scene))
            except Exception as e:
                log.error('[GEOCODE] -- {scene} -- {error}'.format(scene=scene, error=e))
                continue
            
            print('###### SNAP NOISE_POWER: {scene}'.format(scene=scene))
            start_time = time.time()
            try:
                noise_power(infile=scene, outdir=config['out_dir'], polarizations=['VV', 'VH'],
                            spacing=geocode_prms['tr'], t_srs=epsg, refarea='gamma0', tmpdir=config['tmp_dir'],
                            demName=geocode_prms['demName'], externalDEMFile=config.get('ext_dem_file'),
                            externalDEMApplyEGM=geocode_prms['externalDEMApplyEGM'],
                            alignToStandardGrid=geocode_prms['alignToStandardGrid'],
                            standardGridOriginX=align_dict['xmax'], standardGridOriginY=align_dict['ymin'],
                            clean_edges=True)
                t = round((time.time() - start_time), 2)
                log.info('[NOISE_P] -- {scene} -- {time}'.format(scene=scene, time=t))
            except Exception as e:
                log.error('[NOISE_P] -- {scene} -- {error}'.format(scene=scene, error=e))
                continue
    
    ####################################################################################################################
    # NRB - final product generation
    
    if nrb_flag:
        # Group consecutive acquisitions
        selection_grouped = groupbyTime(images=selection, function=seconds, time=60)
        
        for tile in tiles:
            outdir = os.path.join(config['out_dir'], tile)
            os.makedirs(outdir, exist_ok=True)
            
            for scenes in selection_grouped:
                if isinstance(scenes, str):
                    scenes = [scenes]
                print('###### NRB: {tile} | {scenes}'.format(tile=tile, scenes=[os.path.basename(s) for s in scenes]))
                start_time = time.time()
                try:
                    nrb_processing(scenes=scenes, datadir=os.path.dirname(outdir), outdir=outdir, tile=tile,
                                   extent=geo_dict[tile]['ext'],
                                   epsg=epsg, external_wbm=config.get('ext_wbm_file'), dem_name=geocode_prms['demName'])
                    log.info('[    NRB] -- {scenes} -- {time}'.format(scenes=scenes,
                                                                      time=round((time.time() - start_time), 2)))
                except Exception as e:
                    log.error('[    NRB] -- {scenes} -- {error}'.format(scenes=scenes, error=e))
                    continue
