import os
import time
from osgeo import gdal
from spatialist import bbox, intersect
from spatialist.ancillary import finder
from pyroSAR import identify, identify_many, Archive
from s1ard import etad, dem, ard, snap
from s1ard.config import get_config, snap_conf, gdal_conf
import s1ard.ancillary as anc
import s1ard.tile_extraction as tile_ex
from s1ard import search
from s1ard import ocn

gdal.UseExceptions()


def main(config_file, section_name='PROCESSING', debug=False, **kwargs):
    """
    Main function that initiates and controls the processing workflow.
    
    Parameters
    ----------
    config_file: str
        Full path to a `config.ini` file.
    section_name: str
        Section name of the `config.ini` file that processing parameters
        should be parsed from. Default is 'PROCESSING'.
    debug: bool
        Set logging level to DEBUG? Default is False.
    **kwargs
        extra arguments to override parameters in the config file. E.g. `acq_mode`.
    """
    update = False  # update existing products? Internal development flag.
    config = get_config(config_file=config_file, proc_section=section_name, **kwargs)
    log = anc.set_logging(config=config, debug=debug)
    geocode_prms = snap_conf(config=config)
    gdal_prms = gdal_conf(config=config)
    
    anc.check_spacing(geocode_prms['spacing'])
    
    sar_flag = 'sar' in config['mode']
    nrb_flag = 'nrb' in config['mode']
    orb_flag = 'orb' in config['mode']
    
    # DEM download authentication
    username, password = dem.authenticate(dem_type=config['dem_type'],
                                          username=None, password=None)
    ####################################################################################################################
    # scene selection
    log.info('collecting scenes')
    
    if config['db_file'] is not None:
        scenes = finder(config['scene_dir'], [r'^S1[AB].*(SAFE|zip)$'],
                        regex=True, recursive=True, foldermode=1)
        archive = Archive(dbfile=config['db_file'])
        archive.insert(scenes)
    else:
        archive = search.STACArchive(url=config['stac_catalog'],
                                     collections=config['stac_collections'])
    
    if config['scene'] is None:
        attr_search = ['sensor', 'product', 'mindate', 'maxdate',
                       'aoi_tiles', 'aoi_geometry', 'date_strict']
        dict_search = {k: config[k] for k in attr_search}
        dict_search['acquisition_mode'] = config['acq_mode']
        
        if config['datatake'] is not None:
            frame_number = [int(x, 16) for x in config['datatake']]
        else:
            frame_number = None
        dict_search['frameNumber'] = frame_number
        
        selection, aoi_tiles = search.scene_select(archive=archive,
                                                   **dict_search)
        
        if len(selection) == 0:
            log.error('could not find any scenes')
            archive.close()
            return
        
        log.info(f'found {len(selection)} scene(s)')
        scenes = identify_many(selection, sortkey='start')
        search.check_acquisition_completeness(scenes=scenes, archive=archive)
    else:
        if config['mode'] != ['sar']:
            raise RuntimeError("if argument 'scene' is set, the processing mode must be 'sar'")
        scenes = [identify(config['scene'])]
        aoi_tiles = []
    ####################################################################################################################
    # get neighboring GRD scenes to add a buffer to the geocoded scenes
    # otherwise there will be a gap between final geocoded images.
    if config['product'] == 'GRD' and sar_flag:
        log.info('collecting GRD neighbors')
        neighbors = []
        for scene in scenes:
            neighbors.append(search.collect_neighbors(archive=archive, scene=scene))
    else:
        neighbors = [None] * len(scenes)
    ####################################################################################################################
    # OCN scene selection
    if 'wm' in config['annotation']:
        log.info('collecting OCN products')
        scenes_ocn = []
        for scene in scenes:
            start, stop = anc.buffer_time(scene.start, scene.stop, seconds=2)
            result = archive.select(product='OCN', mindate=start,
                                    maxdate=stop, date_strict=True)
            if len(result) == 1:
                scenes_ocn.append(result[0])
            else:
                if len(result) == 0:
                    log.error(f'could not find an OCN product for scene {scene.scene}')
                else:
                    msg = f'found multiple OCN products for scene {scene.scene}:\n'
                    log.error(msg + '\n'.join(result))
                archive.close()
                return
        scenes_ocn = identify_many(scenes_ocn)
    else:
        scenes_ocn = []
    
    archive.close()
    ####################################################################################################################
    # annotation layer selection
    annotation = config.get('annotation', None)
    measurement = config['measurement']
    export_extra = None
    lookup = {'dm': 'layoverShadowMask',
              'ei': 'incidenceAngleFromEllipsoid',
              'lc': 'scatteringArea',
              'ld': 'lookDirection',
              'li': 'localIncidenceAngle',
              'np': 'NESZ',
              'gs': 'gammaSigmaRatio',
              'sg': 'sigmaGammaRatio'}
    
    if annotation is not None:
        annotation = ['gs' if x == 'ratio' and measurement == 'gamma' else 'sg' if x == 'ratio'
        else x for x in annotation]
        export_extra = []
        for layer in annotation:
            if layer in lookup:
                export_extra.append(lookup[layer])
    ####################################################################################################################
    # main SAR processing
    if sar_flag:
        for i, scene in enumerate(scenes):
            scene_base = os.path.splitext(os.path.basename(scene.scene))[0]
            out_dir_scene = os.path.join(config['sar_dir'], scene_base)
            tmp_dir_scene = os.path.join(config['tmp_dir'], scene_base)
            
            log.info(f'SAR processing scene {i + 1}/{len(scenes)}: {scene.scene}')
            if os.path.isdir(out_dir_scene) and not update:
                log.info('Already processed - Skip!')
                continue
            else:
                os.makedirs(out_dir_scene, exist_ok=True)
                os.makedirs(tmp_dir_scene, exist_ok=True)
            ############################################################################################################
            # Preparation of DEM for SAR processing
            dem_type_lookup = {'Copernicus 10m EEA DEM': 'EEA10',
                               'Copernicus 30m Global DEM II': 'GLO30II',
                               'Copernicus 30m Global DEM': 'GLO30',
                               'GETASSE30': 'GETASSE30'}
            dem_type_short = dem_type_lookup[config['dem_type']]
            fname_base_dem = scene_base + f'_DEM_{dem_type_short}.tif'
            fname_dem = os.path.join(tmp_dir_scene, fname_base_dem)
            os.makedirs(tmp_dir_scene, exist_ok=True)
            if not os.path.isfile(fname_dem):
                log.info('creating scene-specific DEM mosaic')
                with scene.bbox() as geom:
                    dem.mosaic(geometry=geom, outname=fname_dem, dem_type=config['dem_type'],
                               username=username, password=password)
            else:
                log.info(f'found scene-specific DEM mosaic: {fname_dem}')
            ############################################################################################################
            # ETAD correction
            if config['etad']:
                log.info('ETAD correction')
                scene = etad.process(scene=scene, etad_dir=config['etad_dir'],
                                     out_dir=tmp_dir_scene)
            ############################################################################################################
            # determination of look factors
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
            ############################################################################################################
            # main processing routine
            start_time = time.time()
            try:
                log.info('starting SNAP processing')
                snap.process(scene=scene.scene, outdir=config['sar_dir'],
                             measurement=measurement,
                             tmpdir=config['tmp_dir'],
                             dem=fname_dem, neighbors=neighbors[i],
                             export_extra=export_extra,
                             gpt_args=config['snap_gpt_args'],
                             rlks=rlks, azlks=azlks, **geocode_prms)
                t = round((time.time() - start_time), 2)
                log.info(f'SNAP processing finished in {t} seconds')
            except Exception as e:
                log.error(msg=e)
                raise
    ####################################################################################################################
    # OCN preparation
    if len(scenes_ocn) > 0:
        log.info('extracting OCN products')
        for scene in scenes_ocn:
            if scene.compression is not None:
                scene.unpack(directory=config['tmp_dir'], exist_ok=True)
            basename = os.path.basename(scene.scene).replace('.SAFE', '')
            outdir = os.path.join(config['sar_dir'], basename)
            os.makedirs(outdir, exist_ok=True)
            for v in ['owiNrcsCmod', 'owiEcmwfWindSpeed', 'owiEcmwfWindDirection']:
                out = os.path.join(outdir, f'{v}.tif')
                if not os.path.isfile(out):
                    ocn.extract(src=scene.scene, dst=out, variable=v)
    ####################################################################################################################
    # ARD - final product generation
    if nrb_flag or orb_flag:
        product_type = 'NRB' if nrb_flag else 'ORB'
        
        # prepare WBM MGRS tiles
        log.info('preparing WBM tiles')
        vec = [x.geometry() for x in scenes]
        extent = anc.get_max_ext(geometries=vec)
        del vec
        with bbox(coordinates=extent, crs=4326) as box:
            dem.prepare(vector=box, threads=gdal_prms['threads'],
                        dem_dir=None, wbm_dir=config['wbm_dir'],
                        dem_type=config['dem_type'],
                        tilenames=aoi_tiles, username=username, password=password,
                        dem_strict=True)
        log.info(f'grouping scenes by time and starting {product_type} production')
        selection_grouped = anc.group_by_time(scenes=scenes)
        for s, group in enumerate(selection_grouped):
            gstr = f'{s + 1}/{len(selection_grouped)}'
            # check that the scenes can really be grouped together
            anc.check_scene_consistency(scenes=group)
            # get the geometries of all tiles that overlap with the current scene group
            vec = [x.geometry() for x in group]
            tiles = tile_ex.tile_from_aoi(vector=vec,
                                          return_geometries=True,
                                          tilenames=aoi_tiles)
            del vec
            t_total = len(tiles)
            for t, tile in enumerate(tiles):
                # select all scenes from the group whose footprint overlaps with the current tile
                scenes_sub = [x for x in group if intersect(tile, x.geometry())]
                scenes_sub_fnames = [x.scene for x in scenes_sub]
                outdir = os.path.join(config['ard_dir'], tile.mgrs)
                os.makedirs(outdir, exist_ok=True)
                fname_wbm = os.path.join(config['wbm_dir'], config['dem_type'],
                                         '{}_WBM.tif'.format(tile.mgrs))
                if not os.path.isfile(fname_wbm):
                    fname_wbm = None
                add_dem = True  # add the DEM as output layer?
                dem_type = config['dem_type'] if add_dem else None
                extent = tile.extent
                epsg = tile.getProjection('epsg')
                log.info(f'creating product {t + 1}/{t_total} of group {gstr}')
                log.info(f'selected scene(s): {scenes_sub_fnames}')
                try:
                    msg = ard.format(config=config, product_type=product_type,
                                     scenes=scenes_sub_fnames, datadir=config['sar_dir'],
                                     outdir=outdir, tile=tile.mgrs, extent=extent, epsg=epsg,
                                     wbm=fname_wbm, dem_type=dem_type,
                                     multithread=gdal_prms['multithread'], annotation=annotation,
                                     update=update)
                    if msg == 'Already processed - Skip!':
                        log.info(msg)
                except Exception as e:
                    log.error(msg=e)
                    raise
            del tiles
        gdal.SetConfigOption('GDAL_NUM_THREADS', gdal_prms['threads_before'])
