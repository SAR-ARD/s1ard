import os
import re
import time
import numpy as np
from osgeo import gdal
from spatialist import Raster
from spatialist.ancillary import finder
from spatialist.auxil import gdalwarp
from pyroSAR import identify, Archive
from pyroSAR.snap import geocode
from pyroSAR.ancillary import groupbyTime, seconds, find_datasets

from S1_NRB.config import get_config, geocode_params
import S1_NRB.ancillary as ancil
import S1_NRB.tile_extraction as tile_ex
from S1_NRB.metadata import extract, xmlparser, stacparser

gdal.UseExceptions()


def nrb_processing(scenes, outdir, tile, extent, epsg, dem_name, compress='LERC_ZSTD', overviews=None, recursive=False,
                   external_wbm=None):
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
    
    id = identify(scenes[0])
    
    start = id.start
    if len(scenes) == 2:
        id2 = identify(scenes[1])
        start2 = id2.start
        stop = id2.stop
    else:
        start2 = None
        stop = id.stop
    
    meta = {'mission': id.sensor,
            'mode': id.meta['acquisition_mode'],
            'start': start,
            'orbitnumber': id.meta['orbitNumbers_abs']['start'],
            'datatake': hex(id.meta['frameNumber']).replace('x', '').upper(),
            'stop': stop,
            'tile': tile,
            'id': 'ABCD'}
    
    skeleton = '{mission}_{mode}_NRB__1SDV_{start}_{stop}_{orbitnumber:06}_{datatake}_{tile}_{id}'
    
    nrbdir = os.path.join(outdir, skeleton.format(**meta))
    os.makedirs(nrbdir, exist_ok=True)
    
    metaL = meta.copy()
    for key, val in metaL.items():
        if not isinstance(val, int):
            metaL[key] = val.lower()
    
    files1 = find_datasets(directory=os.path.dirname(outdir), recursive=recursive,
                           start=start, stop=start)
    if start2 is not None:
        files2 = find_datasets(directory=os.path.dirname(outdir), recursive=recursive,
                               start=start2, stop=start2)
        files = list(zip(files1, files2))
    else:
        files = files1
    
    if len(files) == 0:
        raise RuntimeError("No pyroSAR datasets were found in the directory '{}'".format(os.path.dirname(outdir)))
    
    skeleton = '{mission}-{mode}-nrb-{start}-{stop}-{orbitnumber:06}-{datatake}-{tile}-{suffix}.tif'
    
    suffices = {'VV': 'vv-g-lin',
                'VH': 'vh-g-lin',
                'HH': 'hh-g-lin',
                'HV': 'hv-g-lin',
                'incidenceAngleFromEllipsoid': 'ei',
                'layoverShadowMask': 'dm',
                'localIncidenceAngle': 'li',
                'scatteringArea': 'lc',
                'gammaSigmaRatio': 'gs'}
    
    # Maximum error threshold on values for LERC* compression.
    # Will be ignored if a compression algorithm is used that isn't related to LERC.
    z_errors = {'VV': 1e-4,
                'VH': 1e-4,
                'HH': 1e-4,
                'HV': 1e-4,
                'incidenceAngleFromEllipsoid': 1e-3,
                'layoverShadowMask': 0,
                'localIncidenceAngle': 1e-2,
                'scatteringArea': 0.1,
                'noisePower': 2e-5,
                'gammaSigmaRatio': 1e-4}
    
    driver = 'COG'
    write_options_base = ['BLOCKSIZE=512', 'OVERVIEW_RESAMPLING=AVERAGE']
    write_options = dict()
    for key in z_errors:
        write_options[key] = write_options_base.copy()
        if compress is not None:
            entry = 'COMPRESS={}'.format(compress)
            write_options[key].append(entry)
            if compress.startswith('LERC'):
                entry = 'MAX_Z_ERROR={:f}'.format(z_errors[key])
                write_options[key].append(entry)

    ####################################################################################################################
    # format existing datasets found by `pyroSAR.ancillary.find_datasets`
    
    pattern = '|'.join(suffices.keys())
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
        
        val = suffices[key]
        metaL['suffix'] = val
        outname_base = skeleton.format(**metaL)
        if re.search('[HV]{2}', key):
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
    ancil.modify_data_mask(dm_path=dm_path, extent=extent, epsg=epsg, driver=driver,
                           creation_opt=write_options['layoverShadowMask'], overviews=overviews,
                           multilayer=True, wbm=wbm, wbm_path=external_wbm)
    
    ####################################################################################################################
    # noise power images
    
    for item in measure:
        pol = re.search('[hv]{2}', item).group()
        np_name = item.replace('measurement', 'annotation').replace(pol + '-g-lin', 'np-{}'.format(pol))
        
        if not os.path.isfile(np_name):
            print(np_name)
            with Raster(item) as ras:
                out = np.random.rand(ras.cols, ras.rows)
                with Raster(measure[0]) as ref:
                    ref_arr = ref.array()
                out[np.isnan(ref_arr)] = np.nan
                ras.write(np_name, format=driver, array=out.astype('float32'),
                          overviews=overviews, options=write_options['noisePower'])
    
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
    
    meta = extract.meta_dict(target=nrbdir, sources=scenes, dem_name=dem_name, proc_time=proc_time)
    xmlparser.main(meta=meta, target=nrbdir, sources=scenes)
    stacparser.main(meta=meta, target=nrbdir)


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
    
    scenes = finder(config['work_dir'], [r'^S1[AB].*\.zip'], regex=True, recursive=True)
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
    # geocode - SNAP processing
    
    if geocode_flag:
        # Process scenes individually
        for scene in selection:
            print('###### SNAP GEOCODE: {scene}'.format(scene=scene))
            start_time = time.time()
            try:
                geocode(infile=scene, outdir=config['out_dir'], t_srs=epsg, tmpdir=config['tmp_dir'],
                        standardGridOriginX=align_dict['xmax'], standardGridOriginY=align_dict['ymin'],
                        externalDEMFile=config['ext_dem_file'], **geocode_prms)
                
                t = round((time.time() - start_time), 2)
                log.info('[GEOCODE] -- {scene} -- {time}'.format(scene=scene, time=t))
                if t <= 500:
                    log.warning('[GEOCODE] -- {scene} -- Processing might have terminated prematurely. Check'
                                ' terminal for uncaught SNAP errors!'.format(scene=scene))
            except Exception as e:
                log.error('[GEOCODE] -- {scene} -- {error}'.format(scene=scene, error=e))
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
                    nrb_processing(scenes=scenes, outdir=outdir, tile=tile, extent=geo_dict[tile]['ext'],
                                   epsg=epsg, external_wbm=config['ext_wbm_file'], dem_name=geocode_prms['demName'])
                    log.info('[    NRB] -- {scenes} -- {time}'.format(scenes=scenes,
                                                                      time=round((time.time() - start_time), 2)))
                except Exception as e:
                    log.error('[    NRB] -- {scenes} -- {error}'.format(scenes=scenes, error=e))
                    continue
