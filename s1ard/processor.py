import os
import time
import inspect
from osgeo import gdal
from spatialist import bbox, intersect
from spatialist.ancillary import finder
from pyroSAR import identify, identify_many, Archive
from s1ard import etad, ard
from s1ard.config import get_config, gdal_conf
from s1ard.ancillary import set_logging
from s1ard import search
from s1ard import ocn
from cesard import dem
import cesard.tile_extraction as tile_ex
from cesard.search import scene_select
from cesard.ancillary import (buffer_time, check_scene_consistency,
                              check_spacing, get_max_ext, group_by_attr)

from s1ard.processors.registry import load_processor

gdal.UseExceptions()


def main(config_file=None, debug=False, **kwargs):
    """
    Main function that initiates and controls the processing workflow.
    
    Parameters
    ----------
    config_file: str or None
        Full path to a `config.ini` file or `None` to use the package's default file.
    debug: bool
        Set logging level to DEBUG? Default is False.
    **kwargs
        extra arguments to override parameters in the config file. E.g. `acq_mode`.
    """
    update = False  # update existing products? Internal development flag.
    config = get_config(config_file=config_file, **kwargs)
    log = set_logging(config=config, debug=debug)
    config_proc = config['processing']
    processor_name = config_proc['processor']
    processor = load_processor(processor_name)
    config_sar = config[processor_name]
    gdal_prms = gdal_conf(config=config)
    
    spacings = {'IW': 10, 'SM': 10, 'EW': 30}
    config_sar['spacing'] = spacings[config_proc['acq_mode']]
    
    check_spacing(config_sar['spacing'])
    
    sar_flag = 'sar' in config_proc['mode']
    nrb_flag = 'nrb' in config_proc['mode']
    orb_flag = 'orb' in config_proc['mode']
    
    # DEM download authentication
    username, password = dem.authenticate(dem_type=config_proc['dem_type'],
                                          username=None, password=None)
    ####################################################################################################################
    # scene selection
    log.info('collecting scenes')
    
    db_file_set = config_proc['db_file'] is not None
    scene_dir_set = config_proc['scene_dir'] is not None
    stac_catalog_set = config_proc['stac_catalog'] is not None
    stac_collections_set = config_proc['stac_collections'] is not None
    parquet_set = config_proc['parquet'] is not None
    
    if db_file_set:
        archive = Archive(dbfile=config_proc['db_file'])
        if scene_dir_set:
            scenes = finder(target=config_proc['scene_dir'],
                            matchlist=[r'^S1[ABCD].*(SAFE|zip)$'],
                            regex=True, recursive=True, foldermode=1)
            archive.insert(scenes)
    elif stac_catalog_set and stac_collections_set:
        archive = search.STACArchive(url=config_proc['stac_catalog'],
                                     collections=config_proc['stac_collections'])
    elif parquet_set:
        archive = search.STACParquetArchive(files=config_proc['parquet'])
    else:
        raise RuntimeError('could not select a search option. Please check your configuration.')
    
    if config_proc['scene'] is None:
        attr_search = ['sensor', 'product', 'mindate', 'maxdate',
                       'aoi_tiles', 'aoi_geometry', 'date_strict']
        dict_search = {k: config_proc[k] for k in attr_search}
        dict_search['acquisition_mode'] = config_proc['acq_mode']
        
        if config_proc['datatake'] is not None:
            frame_number = [int(x, 16) for x in config_proc['datatake']]
        else:
            frame_number = None
        dict_search['frameNumber'] = frame_number
        
        selection, aoi_tiles = scene_select(archive=archive,
                                            **dict_search)
        
        if len(selection) == 0:
            log.error('could not find any scenes')
            archive.close()
            return
        
        log.info(f'found {len(selection)} scene(s)')
        scenes = identify_many(selection, sortkey='start')
        search.check_acquisition_completeness(scenes=scenes, archive=archive)
    else:
        if config_proc['mode'] != ['sar']:
            raise RuntimeError("if argument 'scene' is set, the processing mode must be 'sar'")
        scenes = [identify(config_proc['scene'])]
        config_proc['acq_mode'] = scenes[0].acquisition_mode
        config_proc['product'] = scenes[0].product
        aoi_tiles = []
    
    # group scenes by datatake
    scenes_grouped = group_by_attr(scenes, lambda x: x.meta['frameNumber'])
    
    for scenes in scenes_grouped:
        # check that the scenes can really be grouped together
        check_scene_consistency(scenes=scenes)
        
        # Remove scenes with an invalid (0) slice number if others have a valid one (>0).
        # This ensures that scenes with a valid slice number are preferred.
        slice_numbers = [x.meta['sliceNumber'] for x in scenes]
        if None not in slice_numbers:
            if min(slice_numbers) == 0 and max(slice_numbers) > 0:
                for i in reversed(range(len(scenes))):
                    if slice_numbers[i] == 0:
                        del scenes[i]
                        del slice_numbers[i]
                for i in range(1, len(scenes)):
                    if slice_numbers[i] != slice_numbers[i - 1] + 1:
                        raise RuntimeError(f"nonconsecutive scene group, "
                                           f"slice numbers: {slice_numbers}")
    ####################################################################################################################
    # Get neighboring GRD scenes to add a buffer to the geocoded scenes.
    # Otherwise, there will be a gap between final geocoded images.
    # Buffering is only possible if the product composition is 'Sliced'
    # (not 'Assembled' or 'Individual') and thus has a sliceNumber attribute.
    if config_proc['product'] == 'GRD' and sar_flag:
        log.info('collecting GRD neighbors')
        neighbors = []
        for scenes in scenes_grouped:
            if scenes[0].meta['sliceNumber'] is not None:
                neighbors_group = []
                for scene in scenes:
                    neighbors_scene = search.collect_neighbors(archive=archive,
                                                               scene=scene)
                    neighbors_group.append(neighbors_scene)
            else:
                neighbors_group = [None] * len(scenes)
            neighbors.append(neighbors_group)
    else:
        neighbors = [[None] * len(x) for x in scenes_grouped]
    ####################################################################################################################
    # OCN scene selection
    if 'wm' in config_proc['annotation']:
        log.info('collecting OCN products')
        scenes_ocn = []
        for scenes in scenes_grouped:
            scenes_ocn_group = []
            for scene in scenes:
                start, stop = buffer_time(scene.start, scene.stop, seconds=2)
                result = archive.select(product='OCN', mindate=start,
                                        maxdate=stop, date_strict=True)
                if len(result) == 1:
                    scenes_ocn_group.append(result[0])
                else:
                    if len(result) == 0:
                        log.error(f'could not find an OCN product for scene {scene.scene}')
                    else:
                        msg = f'found multiple OCN products for scene {scene.scene}:\n'
                        log.error(msg + '\n'.join(result))
                    archive.close()
                    return
            scenes_ocn_group = identify_many(scenes_ocn_group)
            scenes_ocn.append(scenes_ocn_group)
    else:
        scenes_ocn = [[] for scenes in scenes_grouped]
    
    archive.close()
    ####################################################################################################################
    # annotation layer selection
    annotation = config_proc['annotation']
    measurement = config_proc['measurement']
    export_extra = processor.translate_annotation(annotation=annotation,
                                                  measurement=measurement)
    ####################################################################################################################
    # main SAR processing
    if sar_flag:
        for h, scenes in enumerate(scenes_grouped):
            log.info(f'SAR processing of group {h + 1}/{len(scenes_grouped)}')
            for i, scene in enumerate(scenes):
                scene_base = os.path.splitext(os.path.basename(scene.scene))[0]
                out_dir_scene = os.path.join(config_proc['sar_dir'], scene_base)
                tmp_dir_scene = os.path.join(config_proc['tmp_dir'], scene_base)
                
                log.info(f'processing scene {i + 1}/{len(scenes)}: {scene.scene}')
                if os.path.isdir(out_dir_scene) and not update:
                    log.info('Already processed - Skip!')
                    continue
                else:
                    os.makedirs(out_dir_scene, exist_ok=True)
                    os.makedirs(tmp_dir_scene, exist_ok=True)
                ########################################################################################################
                # Preparation of DEM for SAR processing
                dem_prepare_mode = config_sar['dem_prepare_mode']
                if dem_prepare_mode is not None:
                    fname_dem = dem.prepare(scene=scene, dem_type=config_proc['dem_type'],
                                            dir_out=tmp_dir_scene, username=username,
                                            password=password, mode=dem_prepare_mode,
                                            tr=(config_sar['spacing'], config_sar['spacing']))
                else:
                    fname_dem = None
                ########################################################################################################
                # ETAD correction
                if config_proc['etad']:
                    log.info('ETAD correction')
                    scene = etad.process(scene=scene, etad_dir=config_proc['etad_dir'],
                                         out_dir=tmp_dir_scene)
                ########################################################################################################
                # determination of look factors
                if scene.product == 'SLC':
                    rlks = {'IW': 5,
                            'SM': 6,
                            'EW': 3}[config_proc['acq_mode']]
                    rlks *= int(config_sar['spacing'] / 10)
                    azlks = {'IW': 1,
                             'SM': 6,
                             'EW': 1}[config_proc['acq_mode']]
                    azlks *= int(config_sar['spacing'] / 10)
                else:
                    rlks = azlks = None
                ########################################################################################################
                # main processing routine
                start_time = time.time()
                try:
                    log.info('starting SAR processing')
                    proc_args = {'scene': scene.scene,
                                 'outdir': config_proc['sar_dir'],
                                 'measurement': measurement,
                                 'tmpdir': config_proc['tmp_dir'],
                                 'dem': fname_dem,
                                 'neighbors': neighbors[h][i],
                                 'export_extra': export_extra,
                                 'rlks': rlks, 'azlks': azlks}
                    proc_args.update(config_sar)
                    sig = inspect.signature(processor.process)
                    accepted_params = set(sig.parameters.keys())
                    proc_args = {k: v for k, v in proc_args.items() if k in accepted_params}
                    processor.process(**proc_args)
                    t = round((time.time() - start_time), 2)
                    log.info(f'SAR processing finished in {t} seconds')
                except Exception as e:
                    log.error(msg=e)
                    raise
    ####################################################################################################################
    # OCN preparation
    if sum(len(x) for x in scenes_ocn) > 0:
        log.info('extracting OCN products')
        for scenes in scenes_ocn:
            for scene in scenes:
                if scene.compression is not None:
                    scene.unpack(directory=config_proc['tmp_dir'], exist_ok=True)
                basename = os.path.basename(scene.scene).replace('.SAFE', '')
                outdir = os.path.join(config_proc['sar_dir'], basename)
                os.makedirs(outdir, exist_ok=True)
                for v in ['owiNrcsCmod', 'owiEcmwfWindSpeed', 'owiEcmwfWindDirection']:
                    out = os.path.join(outdir, f'{v}.tif')
                    if not os.path.isfile(out):
                        ocn.extract(src=scene.scene, dst=out, variable=v)
    ####################################################################################################################
    # ARD - final product generation
    if nrb_flag or orb_flag:
        product_type = 'NRB' if nrb_flag else 'ORB'
        log.info(f'starting {product_type} production')
        
        for s, scenes in enumerate(scenes_grouped):
            log.info(f'ARD processing of group {s + 1}/{len(scenes_grouped)}')
            log.info('preparing WBM tiles')
            vec = [x.geometry() for x in scenes]
            extent = get_max_ext(geometries=vec)
            with bbox(coordinates=extent, crs=4326) as box:
                dem.retile(vector=box, threads=gdal_prms['threads'],
                           dem_dir=None, wbm_dir=config_proc['wbm_dir'],
                           dem_type=config_proc['dem_type'],
                           tilenames=aoi_tiles, username=username, password=password,
                           dem_strict=True)
            # get the geometries of all tiles that overlap with the current scene group
            tiles = tile_ex.tile_from_aoi(vector=vec,
                                          return_geometries=True,
                                          tilenames=aoi_tiles)
            del vec
            t_total = len(tiles)
            for t, tile in enumerate(tiles):
                # select all scenes from the group whose footprint overlaps with the current tile
                scenes_sub = [x for x in scenes if intersect(tile, x.geometry())]
                scenes_sub_fnames = [x.scene for x in scenes_sub]
                fname_wbm = os.path.join(config_proc['wbm_dir'],
                                         config_proc['dem_type'],
                                         '{}_WBM.tif'.format(tile.mgrs))
                if not os.path.isfile(fname_wbm):
                    fname_wbm = None
                add_dem = True  # add the DEM as output layer?
                dem_type = config_proc['dem_type'] if add_dem else None
                extent = tile.extent
                epsg = tile.getProjection('epsg')
                log.info(f'creating product {t + 1}/{t_total}')
                log.info(f'selected scene(s): {scenes_sub_fnames}')
                try:
                    prod_meta = ard.product_info(product_type=product_type, src_ids=scenes_sub,
                                                 tile_id=tile.mgrs, extent=extent, epsg=epsg,
                                                 dir_ard=config_proc['ard_dir'], update=update)
                except RuntimeError:
                    log.info('Already processed - Skip!')
                    del tiles
                    return
                log.info(f'product name: {prod_meta['dir_ard_product']}')
                try:
                    src_ids, sar_assets = ard.get_datasets(scenes=scenes_sub_fnames,
                                                           sar_dir=config_proc['sar_dir'],
                                                           extent=extent, epsg=epsg,
                                                           processor_name=processor_name)
                    
                    ard_assets = ard.format(config=config, prod_meta=prod_meta,
                                            src_ids=src_ids, sar_assets=sar_assets,
                                            tile=tile.mgrs, extent=extent, epsg=epsg,
                                            wbm=fname_wbm, dem_type=dem_type, compress='LERC_ZSTD',
                                            multithread=gdal_prms['multithread'], annotation=annotation)
                    
                    if ard_assets is not None:
                        ard.append_metadata(config=config, prod_meta=prod_meta,
                                            src_ids=src_ids, assets=ard_assets, compression='LERC_ZSTD')
                
                except Exception as e:
                    log.error(msg=e)
                    raise
            del tiles
        gdal.SetConfigOption('GDAL_NUM_THREADS', gdal_prms['threads_before'])
