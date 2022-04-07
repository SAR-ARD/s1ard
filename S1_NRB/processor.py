import os
import time
from osgeo import gdal
from spatialist.ancillary import finder
from pyroSAR import identify_many, Archive
from pyroSAR.snap.util import geocode, noise_power
from pyroSAR.ancillary import groupbyTime, seconds
from S1_NRB import etad, dem, nrb
from S1_NRB.config import get_config, geocode_conf, gdal_conf
import S1_NRB.ancillary as ancil
import S1_NRB.tile_extraction as tile_ex

gdal.UseExceptions()


def main(config_file, section_name, debug=False):
    config = get_config(config_file=config_file, section_name=section_name)
    log = ancil.set_logging(config=config, debug=debug)
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
    
    with Archive(dbfile=config['db_file']) as archive:
        archive.insert(scenes)
        selection = archive.select(product='SLC',
                                   acquisition_mode=config['acq_mode'],
                                   mindate=config['mindate'], maxdate=config['maxdate'])
    ids = identify_many(selection)
    
    if len(selection) == 0:
        message = "No scenes could be found for acquisition mode '{acq_mode}', " \
                  "mindate '{mindate}' and maxdate '{maxdate}' in directory '{scene_dir}'."
        raise RuntimeError(message.format(acq_mode=config['acq_mode'], mindate=config['mindate'],
                                          maxdate=config['maxdate'], scene_dir=config['scene_dir']))
    ####################################################################################################################
    # general setup
    geo_dict = tile_ex.get_tile_dict(config=config, spacing=geocode_prms['spacing'])
    aoi_tiles = list(geo_dict.keys())
    aoi_tiles.remove('align')
    
    epsg_set = set([geo_dict[tile]['epsg'] for tile in list(geo_dict.keys()) if tile != 'align'])
    if len(epsg_set) != 1:
        raise RuntimeError('The AOI covers multiple UTM zones: {}\n '
                           'This is currently not supported. '
                           'Please refine your AOI.'.format(list(epsg_set)))
    epsg = epsg_set.pop()
    
    np_dict = {'sigma0': 'NESZ', 'beta0': 'NEBZ', 'gamma0': 'NEGZ'}
    np_refarea = 'sigma0'
    ####################################################################################################################
    # DEM download and MGRS-tiling
    geometries = [scene.bbox() for scene in ids]
    dem.prepare(geometries=geometries, threads=gdal_prms['threads'],
                epsg=epsg, spacing=geocode_prms['spacing'],
                dem_dir=config['dem_dir'], wbm_dir=config['wbm_dir'],
                dem_type=config['dem_type'], kml_file=config['kml_file'])
    del geometries
    ####################################################################################################################
    # SNAP RTC processing
    if snap_flag:
        for i, scene in enumerate(ids):
            ###############################################
            scene_base = os.path.splitext(os.path.basename(scene.scene))[0]
            out_dir_scene = os.path.join(config['rtc_dir'], scene_base)
            tmp_dir_scene = os.path.join(config['tmp_dir'], scene_base)
            out_dir_scene_epsg = os.path.join(out_dir_scene, str(epsg))
            tmp_dir_scene_epsg = os.path.join(tmp_dir_scene, str(epsg))
            os.makedirs(out_dir_scene_epsg, exist_ok=True)
            os.makedirs(tmp_dir_scene_epsg, exist_ok=True)
            fname_dem = os.path.join(tmp_dir_scene_epsg,
                                     scene.outname_base() + '_DEM_{}.tif'.format(epsg))
            ###############################################
            # scene-specific DEM preparation
            with scene.bbox() as geometry:
                dem.mosaic(geometry, outname=fname_dem, epsg=epsg,
                           dem_type=config['dem_type'], kml_file=config['kml_file'],
                           dem_dir=config['dem_dir'])
            
            if config['dem_type'] == 'Copernicus 30m Global DEM':
                ex_dem_nodata = -99
            else:
                ex_dem_nodata = None
            ###############################################
            # ETAD correction
            if config['etad']:
                print('###### [   ETAD] Scene {s}/{s_total}: {scene}'.format(s=i + 1, s_total=len(ids),
                                                                             scene=scene.scene))
                scene = etad.process(scene=scene, etad_dir=config['etad_dir'],
                                     out_dir=tmp_dir_scene, log=log)
            ###############################################
            list_processed = finder(out_dir_scene_epsg, ['*'])
            exclude = list(np_dict.values())
            print('###### [GEOCODE] Scene {s}/{s_total}: {scene}'.format(s=i + 1, s_total=len(ids),
                                                                         scene=scene.scene))
            if len([item for item in list_processed if not any(ex in item for ex in exclude)]) < 4:
                start_time = time.time()
                try:
                    geocode(infile=scene, outdir=out_dir_scene_epsg, t_srs=epsg, tmpdir=tmp_dir_scene_epsg,
                            standardGridOriginX=geo_dict['align']['xmax'],
                            standardGridOriginY=geo_dict['align']['ymin'],
                            externalDEMFile=fname_dem, externalDEMNoDataValue=ex_dem_nodata, **geocode_prms)
                    t = round((time.time() - start_time), 2)
                    _log(handler=log, mode='info', proc_step='GEOCODE', scenes=scene.scene, epsg=epsg, msg=t)
                    if t <= 500:
                        msg = 'Processing might have terminated prematurely. Check terminal for uncaught SNAP errors!'
                        _log(handler=log, mode='warning', proc_step='GEOCODE', scenes=scene.scene, epsg=epsg, msg=msg)
                except Exception as e:
                    _log(handler=log, mode='exception', proc_step='GEOCODE', scenes=scene.scene, epsg=epsg, msg=e)
                    continue
            else:
                msg = 'Already processed - Skip!'
                print('### ' + msg)
                _log(handler=log, mode='info', proc_step='GEOCODE', scenes=scene.scene, epsg=epsg, msg=msg)
            ###############################################
            print('###### [NOISE_P] Scene {s}/{s_total}: {scene}'.format(s=i + 1, s_total=len(ids),
                                                                         scene=scene.scene))
            if len([item for item in list_processed if np_dict[np_refarea] in item]) == 0:
                start_time = time.time()
                try:
                    noise_power(infile=scene.scene, outdir=out_dir_scene_epsg, polarizations=scene.polarizations,
                                spacing=geocode_prms['spacing'], refarea=np_refarea, tmpdir=tmp_dir_scene_epsg,
                                externalDEMFile=fname_dem, externalDEMNoDataValue=ex_dem_nodata, t_srs=epsg,
                                externalDEMApplyEGM=geocode_prms['externalDEMApplyEGM'],
                                alignToStandardGrid=geocode_prms['alignToStandardGrid'],
                                standardGridOriginX=geo_dict['align']['xmax'],
                                standardGridOriginY=geo_dict['align']['ymin'],
                                clean_edges=geocode_prms['clean_edges'],
                                clean_edges_npixels=geocode_prms['clean_edges_npixels'],
                                rlks=geocode_prms['rlks'], azlks=geocode_prms['azlks'])
                    _log(handler=log, mode='info', proc_step='NOISE_P', scenes=scene.scene, epsg=epsg,
                         msg=round((time.time() - start_time), 2))
                except Exception as e:
                    _log(handler=log, mode='exception', proc_step='NOISE_P', scenes=scene.scene, epsg=epsg, msg=e)
                    continue
            else:
                msg = 'Already processed - Skip!'
                print('### ' + msg)
                _log(handler=log, mode='info', proc_step='NOISE_P', scenes=scene.scene, epsg=epsg, msg=msg)
    
    ####################################################################################################################
    # NRB - final product generation
    if nrb_flag:
        selection_grouped = groupbyTime(images=selection, function=seconds, time=60)
        for t, tile in enumerate(aoi_tiles):
            outdir = os.path.join(config['nrb_dir'], tile)
            os.makedirs(outdir, exist_ok=True)
            wbm = os.path.join(config['wbm_dir'], config['dem_type'], '{}_WBM.tif'.format(tile))
            if not os.path.isfile(wbm):
                wbm = None
            
            for s, scenes in enumerate(selection_grouped):
                if isinstance(scenes, str):
                    scenes = [scenes]
                print('###### [    NRB] Tile {t}/{t_total}: {tile} | '
                      'Scenes {s}/{s_total}: {scenes} '.format(tile=tile, t=t + 1, t_total=len(aoi_tiles),
                                                               scenes=[os.path.basename(s) for s in scenes],
                                                               s=s + 1, s_total=len(selection_grouped)))
                start_time = time.time()
                try:
                    nrb.format(config=config, scenes=scenes, datadir=config['rtc_dir'], outdir=outdir,
                               tile=tile, extent=geo_dict[tile]['ext'], epsg=epsg, wbm=wbm,
                               multithread=gdal_prms['multithread'])
                    _log(handler=log, mode='info', proc_step='NRB', scenes=scenes, epsg=epsg,
                         msg=round((time.time() - start_time), 2))
                except Exception as e:
                    _log(handler=log, mode='exception', proc_step='NRB', scenes=scenes, epsg=epsg, msg=e)
                    continue
        
        gdal.SetConfigOption('GDAL_NUM_THREADS', gdal_prms['threads_before'])


def _log(handler, mode, proc_step, scenes, epsg, msg):
    """Helper function to format and handle log messages during processing."""
    proc_step = proc_step.zfill(7).replace('0', ' ')
    message = '[{proc_step}] -- {scenes} [{epsg}] -- {msg}'
    log = message.format(proc_step=proc_step, scenes=scenes, epsg=epsg, msg=msg)
    if mode == 'info':
        handler.info(log)
    elif mode == 'warning':
        handler.warning(log)
    elif mode == 'exception':
        handler.exception(log)
    else:
        raise RuntimeError('log mode {} is not supported'.format(mode))
