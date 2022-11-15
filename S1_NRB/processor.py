import os
import time
import itertools
from osgeo import gdal
from spatialist import Vector, bbox
from spatialist.ancillary import finder
from pyroSAR import identify_many, Archive
from S1_NRB import etad, dem, nrb, snap
from S1_NRB.config import get_config, snap_conf, gdal_conf
from S1_NRB.ancillary import set_logging, log, get_max_ext, check_spacing, group_by_time
import S1_NRB.tile_extraction as tile_ex

gdal.UseExceptions()


def main(config_file, section_name='PROCESSING', debug=False):
    """
    Main function that initiates and controls the processing workflow.
    
    Parameters
    ----------
    config_file: str
        Full path to a `config.ini` file.
    section_name: str, optional
        Section name of the `config.ini` file that processing parameters should be parsed from. Default is 'PROCESSING'.
    debug: bool, optional
        Set pyroSAR logging level to DEBUG? Default is False.
    """
    config = get_config(config_file=config_file, proc_section=section_name)
    logger = set_logging(config=config, debug=debug)
    geocode_prms = snap_conf(config=config)
    gdal_prms = gdal_conf(config=config)
    
    check_spacing(geocode_prms['spacing'])
    
    rtc_flag = True
    nrb_flag = True
    if config['mode'] == 'rtc':
        nrb_flag = False
    elif config['mode'] == 'nrb':
        rtc_flag = False
    
    ####################################################################################################################
    # archive / scene selection
    scenes = finder(config['scene_dir'], [r'^S1[AB].*(SAFE|zip)$'],
                    regex=True, recursive=True, foldermode=1)
    
    if config['acq_mode'] == 'SM':
        acq_mode_search = ('S1', 'S2', 'S3', 'S4', 'S5', 'S6')
    else:
        acq_mode_search = config['acq_mode']
    
    vec = [None]
    aoi_tiles = None
    selection = []
    if config['aoi_tiles'] is not None:
        vec = tile_ex.extract_tile(kml=config['kml_file'], tile=config['aoi_tiles'])
        aoi_tiles = config['aoi_tiles']
    elif config['aoi_geometry'] is not None:
        vec = [Vector(config['aoi_geometry'])]
        aoi_tiles = tile_ex.tiles_from_aoi(vectorobject=vec[0], kml=config['kml_file'])
    
    with Archive(dbfile=config['db_file']) as archive:
        archive.insert(scenes)
        for item in vec:
            selection.extend(
                archive.select(vectorobject=item,
                               product=config['product'],
                               acquisition_mode=acq_mode_search,
                               mindate=config['mindate'],
                               maxdate=config['maxdate']))
    selection = list(set(selection))
    del vec
    
    if len(selection) == 0:
        message = "No scenes could be found for acquisition mode '{acq_mode}', " \
                  "mindate '{mindate}' and maxdate '{maxdate}' in directory '{scene_dir}'."
        raise RuntimeError(message.format(acq_mode=config['acq_mode'], mindate=config['mindate'],
                                          maxdate=config['maxdate'], scene_dir=config['scene_dir']))
    scenes = identify_many(selection)
    ####################################################################################################################
    # DEM download and WBM MGRS-tiling
    if rtc_flag:
        username, password = dem.authenticate(dem_type=config['dem_type'], username=None, password=None)
        for i, scene in enumerate(scenes):
            scene_base = os.path.splitext(os.path.basename(scene.scene))[0]
            out_dir_scene = os.path.join(config['rtc_dir'], scene_base)
            tmp_dir_scene = os.path.join(config['tmp_dir'], scene_base)
            
            tiles = tile_ex.tiles_from_aoi(vectorobject=scene.bbox(),
                                           kml=config['kml_file'],
                                           return_geometries=True)
            
            for epsg, group in itertools.groupby(tiles, lambda x: x.getProjection('epsg')):
                geometries = list(group)
                ext = get_max_ext(geometries=geometries, buffer=200)
                with bbox(coordinates=ext, crs=epsg) as geom:
                    geom.reproject(4326)
                    print(f'###### [    DEM] processing EPSG:{epsg}')
                    dem.prepare(geometries=[geom], threads=gdal_prms['threads'],
                                epsg=epsg, dem_dir=None, wbm_dir=config['wbm_dir'],
                                dem_type=config['dem_type'], kml_file=config['kml_file'],
                                username=username, password=password)
            
            dem_type_lookup = {'Copernicus 10m EEA DEM': 'EEA10',
                               'Copernicus 30m Global DEM II': 'GLO30II',
                               'Copernicus 30m Global DEM': 'GLO30',
                               'GETASSE30': 'GETASSE30'}
            dem_type_short = dem_type_lookup[config['dem_type']]
            fname_base_dem = scene_base + f'_DEM_{dem_type_short}.tif'
            fname_dem = os.path.join(tmp_dir_scene, fname_base_dem)
            os.makedirs(tmp_dir_scene, exist_ok=True)
            print('###### [    DEM] creating scene-specific mosaic:', fname_dem)
            with scene.bbox() as geom:
                dem.mosaic(geometry=geom, outname=fname_dem, dem_type=config['dem_type'],
                           username=username, password=password)
            ############################################################################################################
            # RTC processing
            print(f'###### [    RTC] Scene {i + 1}/{len(scenes)}: {scene.scene}')
            if os.path.isdir(out_dir_scene):
                msg = 'Already processed - Skip!'
                print('### ' + msg)
                log(handler=logger, mode='info', proc_step='GEOCODE', scenes=scene.scene, msg=msg)
                continue
            else:
                os.makedirs(out_dir_scene)
                os.makedirs(tmp_dir_scene, exist_ok=True)
            ###############################################
            # ETAD correction
            if config['etad']:
                print('###### [   ETAD] Scene {s}/{s_total}: {scene}'.format(s=i + 1, s_total=len(scenes),
                                                                             scene=scene.scene))
                scene = etad.process(scene=scene, etad_dir=config['etad_dir'],
                                     out_dir=tmp_dir_scene, log=logger)
            ###############################################
            if scene.product == 'SLC':
                rlks = {'IW': 5,
                        'SM': 6,
                        'EW': 3}[config['acq_mode']]
                rlks *= int(geocode_prms['spacing'] / 10)
                azlks = {'IW': 1,
                         'SM': 6,
                         'EW': 1}[config['acq_mode']]
                azlks *= int(geocode_prms['spacing'] / 10)
            else:
                rlks = azlks = None
            ###############################################
            start_time = time.time()
            try:
                snap.process(scene=scene.scene, outdir=config['rtc_dir'],
                             tmpdir=config['tmp_dir'], kml=config['kml_file'],
                             dem=fname_dem,
                             rlks=rlks, azlks=azlks, **geocode_prms)
                t = round((time.time() - start_time), 2)
                log(handler=logger, mode='info', proc_step='RTC', scenes=scene.scene, msg=t)
            except Exception as e:
                log(handler=logger, mode='exception', proc_step='RTC', scenes=scene.scene, msg=e)
                continue
    ####################################################################################################################
    # NRB - final product generation
    if nrb_flag:
        if aoi_tiles is None:
            if config['aoi_tiles'] is not None:
                aoi_tiles = config['aoi_tiles']
            elif config['aoi_geometry'] is not None:
                with Vector(config['aoi_geometry']) as vec:
                    aoi_tiles = tile_ex.tiles_from_aoi(vectorobject=vec, kml=config['kml_file'])
        selection_grouped = group_by_time(scenes=scenes)
        for t, tile in enumerate(aoi_tiles):
            outdir = os.path.join(config['nrb_dir'], tile)
            os.makedirs(outdir, exist_ok=True)
            wbm = os.path.join(config['wbm_dir'], config['dem_type'], '{}_WBM.tif'.format(tile))
            if not os.path.isfile(wbm):
                wbm = None
            
            for s, scenes in enumerate(selection_grouped):
                print('###### [    NRB] Tile {t}/{t_total}: {tile} | '
                      'Scenes {s}/{s_total}: {scenes} '.format(tile=tile, t=t + 1, t_total=len(aoi_tiles),
                                                               scenes=[os.path.basename(s) for s in scenes],
                                                               s=s + 1, s_total=len(selection_grouped)))
                try:
                    with tile_ex.extract_tile(kml=config['kml_file'], tile=tile) as geom:
                        extent = geom.extent
                        epsg = geom.getProjection('epsg')
                    msg = nrb.format(config=config, scenes=scenes, datadir=config['rtc_dir'],
                                     outdir=outdir, tile=tile, extent=extent, epsg=epsg,
                                     wbm=wbm, multithread=gdal_prms['multithread'])
                    if msg == 'Already processed - Skip!':
                        print('### ' + msg)
                    log(handler=logger, mode='info', proc_step='NRB', scenes=scenes, msg=msg)
                except Exception as e:
                    log(handler=logger, mode='exception', proc_step='NRB', scenes=scenes, msg=e)
                    continue
        gdal.SetConfigOption('GDAL_NUM_THREADS', gdal_prms['threads_before'])
