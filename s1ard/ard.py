import os
import re
import shutil
from datetime import datetime, timezone
from copy import deepcopy
from spatialist.vector import Vector, bbox, intersect
from spatialist.raster import Raster
from spatialist.auxil import gdalwarp, gdalbuildvrt
from spatialist.ancillary import finder
from pyroSAR import identify_many
from pyroSAR.ancillary import Lock
from pyroSAR.drivers import ID
import s1ard
from s1ard import ocn
from s1ard.metadata.mapping import ARD_PATTERN
from s1ard.metadata.extract import append_wind_norm, copy_src_meta, meta_dict
from s1ard.processors.registry import load_processor
from cesard.ancillary import (generate_unique_id, vrt_add_overviews,
                              datamask, get_tmp_name)
from cesard.ard import (calc_product_start_stop, create_acq_id_image,
                        create_data_mask, create_rgb_vrt, create_vrt)
from cesard import dem
from cesard.metadata import xml, stac
from cesard.metadata.mapping import LERC_ERR_THRES
import logging

log = logging.getLogger('s1ard')


def product_info(
        product_type: str,
        src_ids: list[ID],
        tile_id: str,
        extent: dict[str, int | float],
        epsg: int,
        dir_ard: str,
        update: bool = False,
        product_id: str | None = None
) -> dict[str, str | int | datetime]:
    """
    Create ARD product metadata.
    
    Parameters
    ----------
    product_type:
        the ARD product type; options: NRB, ORB
    src_ids:
        the source product objects
    tile_id:
        the MGRS tile ID
    extent:
        the extent of the MGRS tile
    epsg:
        the EPSG code of the MGRS tile
    dir_ard:
        the ARD output directory
    update:
        update an existing product (or create a new one)
    product_id:
        an existing product ID. Default None: create a new one using
        function :func:`generate_unique_id`.

    Returns
    -------
        ARD product metadata
    """
    # determine processing timestamp and generate unique ID
    proc_time = datetime.now(timezone.utc)
    if product_id is None:
        t = proc_time.isoformat().encode()
        product_id = generate_unique_id(encoded_str=t)
    
    ard_start, ard_stop = calc_product_start_stop(src_ids=src_ids, extent=extent, epsg=epsg)
    pol_str = '_'.join(sorted(src_ids[0].polarizations))
    meta = {'mission': src_ids[0].sensor,
            'mode': src_ids[0].meta['acquisition_mode'],
            'product_type': product_type,
            'polarization': {'HH': 'SH',
                             'VV': 'SV',
                             'HH_HV': 'DH',
                             'VH_VV': 'DV'}[pol_str],
            'start': ard_start,
            'stop': ard_stop,
            'proc_time': proc_time,
            'orbitnumber': src_ids[0].meta['orbitNumbers_abs']['start'],
            'datatake': hex(src_ids[0].meta['frameNumber']).replace('x', '').upper(),
            'tile': tile_id,
            'id': product_id}
    meta_name = deepcopy(meta)
    meta_name['start'] = datetime.strftime(meta_name['start'], '%Y%m%dT%H%M%S')
    del meta_name['stop']
    del meta_name['proc_time']
    meta_name_lower = dict((k, v.lower() if isinstance(v, str) else v)
                           for k, v in meta_name.items())
    skeleton_dir = ('{mission}_{mode}_{product_type}__1S{polarization}_{start}_'
                    '{orbitnumber:06}_{datatake:0>6}_{tile}_{id}')
    skeleton_files = '{mission}-{mode}-{product_type}-{start}-{orbitnumber:06}-{datatake:0>6}-{tile}'
    
    meta['product_base'] = skeleton_dir.format(**meta_name)
    dir_ard_tile = os.path.join(dir_ard, tile_id)
    os.makedirs(dir_ard_tile, exist_ok=True)
    meta['dir_ard_product'] = os.path.join(dir_ard_tile, meta['product_base'])
    meta['file_base'] = skeleton_files.format(**meta_name_lower) + '-{suffix}.tif'
    
    # check existence of products
    msg = 'Already processed - Skip!'
    pattern = meta['product_base'].replace(product_id, '*')
    existing = finder(dir_ard_tile, [pattern], foldermode=2)
    if len(existing) > 0:
        if not update:
            raise RuntimeError(msg)
        else:
            if existing[0] != meta['dir_ard_product']:
                existing_meta = re.search(ARD_PATTERN, os.path.basename(existing[0])).groupdict()
                return product_info(product_type=product_type, src_ids=src_ids,
                                    tile_id=tile_id, extent=extent, epsg=epsg,
                                    dir_ard=dir_ard, update=update,
                                    product_id=existing_meta['ID'])
            else:
                return meta
    else:
        try:
            os.makedirs(meta['dir_ard_product'], exist_ok=False)
        except OSError:
            raise RuntimeError(msg)
    return meta


def format(
        config: dict[str, dict[str, int | float | str | list[str]]],
        prod_meta: dict[str, str | int | datetime],
        src_ids: list[ID],
        sar_assets: list[dict[str, str]],
        tile: str,
        extent: dict[str, int | float],
        epsg: int,
        wbm: str | None = None,
        dem_type: str | None = None,
        multithread: bool = True,
        compress: str | None = None,
        overviews: list[int] | None = None,
        annotation: list[str] | None = None
) -> list[str] | None:
    """
    Create ARD products from the SAR processor output.
    This includes the following:
    
    - Creating all measurement and annotation datasets in Cloud Optimized GeoTIFF (COG) format
    - Creating additional annotation datasets in Virtual Raster Tile (VRT) format
    - Applying the ARD product directory structure & naming convention
    - Generating metadata in XML and JSON formats for the ARD product as well as source SLC datasets
    
    Parameters
    ----------
    config:
        Dictionary of the parsed config parameters for the current process.
    prod_meta:
        Product metadata as returned by :func:`~s1ard.ard.product_info`.
    src_ids:
        List of scenes to process. Either a single scene or multiple, matching scenes (consecutive acquisitions).
        All scenes are expected to overlap with `extent` and an error will be thrown if the processing output
        cannot be found for any of the scenes.
    sar_assets:
        The SAR processing assets as returned by :func:`get_datasets`.
    tile:
        ID of an MGRS tile.
    extent:
        Spatial extent of the MGRS tile, derived from a :class:`~spatialist.vector.Vector` object.
    epsg:
        The CRS used for the ARD product; provided as an EPSG code.
    wbm:
        Path to a water body mask file with the dimensions of an MGRS tile.
    dem_type:
        if defined, a DEM layer will be added to the product. The suffix `em` (elevation model) is used.
        Default `None`: do not add a DEM layer.
    multithread:
        Should `gdalwarp` use multithreading? Default is True. The number of threads used, can be adjusted in the
        `config.ini` file with the parameter `gdal_threads`.
    compress:
        Compression algorithm to use. See https://gdal.org/drivers/raster/gtiff.html#creation-options for options.
        Defaults to 'LERC_DEFLATE'.
    overviews:
        Internal overview levels to be created for each GeoTIFF file. Defaults to [2, 4, 9, 18, 36]
    annotation:
        an optional list to select the annotation layers. Default `None`: create all layers if the
        source products contain the required input layers. Options:
        
        - dm: data mask (four masks: not layover not shadow, layover, shadow, water)
        - ei: ellipsoidal incident angle
        - em: digital elevation model
        - id: acquisition ID image (source scene ID per pixel)
        - lc: RTC local contributing area
        - ld: range look direction angle
        - li: local incident angle
        - np: noise power (NESZ, per polarization)
        - gs: gamma-sigma ratio: sigma0 RTC / gamma0 RTC
        - sg: sigma-gamma ratio: gamma0 RTC / sigma0 ellipsoidal
        - wm: OCN product wind model; requires OCN scenes via argument `scenes_ocn`
    
    Returns
    -------
        the ARD product assets
    """
    if compress is None:
        compress = 'LERC_ZSTD'
    if overviews is None:
        overviews = [2, 4, 9, 18, 36]
    ovr_resampling = 'AVERAGE'
    driver = 'COG'
    blocksize = 512
    dst_nodata_float = -9999.0
    dst_nodata_byte = 255
    vrt_nodata = 'nan'  # was found necessary for proper calculation of statistics in QGIS
    vrt_options = {'VRTNodata': vrt_nodata}
    
    processor_name = config['processing']['processor']
    
    if len(src_ids) == 0:
        log.error(f'None of the processed scenes overlap with the current tile {tile}')
        shutil.rmtree(prod_meta['dir_ard_product'])
        return
    
    if annotation is not None:
        allowed = []
        for key in sar_assets[0]:
            c1 = re.search('[gs]-lin', key)
            c2 = key in annotation
            c3 = key in ['gs', 'sg'] and 'ratio' in annotation
            c4 = key.startswith('np') and 'np' in annotation
            if c1 or c2 or c3 or c4:
                allowed.append(key)
    else:
        allowed = [key for key in sar_assets[0].keys() if re.search('[gs]-lin', key)]
        annotation = []
    for item in ['em', 'id']:
        if item in annotation:
            allowed.append(item)
    
    # GDAL output bounds
    bounds = [extent['xmin'], extent['ymin'], extent['xmax'], extent['ymax']]
    
    subdirectories = ['measurement', 'annotation', 'source', 'support']
    for subdirectory in subdirectories:
        os.makedirs(os.path.join(prod_meta['dir_ard_product'], subdirectory), exist_ok=True)
    
    # prepare raster write options; https://gdal.org/drivers/raster/cog.html
    write_options_base = ['BLOCKSIZE={}'.format(blocksize),
                          'OVERVIEW_RESAMPLING={}'.format(ovr_resampling)]
    write_options = dict()
    for key in LERC_ERR_THRES:
        write_options[key] = write_options_base.copy()
        if compress is not None:
            entry = 'COMPRESS={}'.format(compress)
            write_options[key].append(entry)
            if compress.startswith('LERC'):
                entry = 'MAX_Z_ERROR={:f}'.format(LERC_ERR_THRES[key])
                write_options[key].append(entry)
    
    # create raster files: linear gamma0/sigma0 backscatter (-[vh|vv|hh|hv]-[gs]-lin.tif),
    # ellipsoidal incident angle (-ei.tif), gamma-to-sigma ratio (-gs.tif),
    # local contributing area (-lc.tif), local incident angle (-li.tif),
    # noise power images (-np-[vh|vv|hh|hv].tif)
    ard_assets = dict()
    for key in list(sar_assets[0].keys()):
        if key in ['dm', 'wm'] or key not in LERC_ERR_THRES.keys() or key not in allowed:
            # raster files for keys 'dm' and 'wm' are created later
            continue
        
        outname_base = prod_meta["file_base"].format(suffix=key)
        if re.search('[gs]-lin', key):
            subdir = 'measurement'
        else:
            subdir = 'annotation'
        outname = os.path.join(prod_meta['dir_ard_product'], subdir, outname_base)
        
        if not os.path.isfile(outname):
            log.info(f"creating {os.path.relpath(outname, prod_meta['dir_ard_product'])}")
            images = [ds[key] for ds in sar_assets]
            ras = None
            if len(images) > 1:
                ras = Raster(images, list_separate=False)
                source = ras.filename
            else:
                source = get_tmp_name(suffix='.vrt')
                gdalbuildvrt(src=images[0], dst=source)
            
            # modify temporary VRT to make sure overview levels and resampling are properly applied
            vrt_add_overviews(vrt=source, overviews=overviews, resampling=ovr_resampling)
            
            options = {'format': driver, 'outputBounds': bounds,
                       'dstNodata': dst_nodata_float, 'multithread': multithread,
                       'creationOptions': write_options[key]}
            
            gdalwarp(src=source, dst=outname, **options)
            if ras is not None:
                ras.close()
        ard_assets[key] = outname
    
    # define a reference raster and list all gamma0/sigma0 backscatter measurement rasters
    measure_tifs = [v for k, v in ard_assets.items() if re.search('[gs]-lin', k)]
    ref_key = list(ard_assets.keys())[0]
    ref_tif = ard_assets[ref_key]
    
    # create data mask raster (-dm.tif)
    if 'dm' in allowed:
        if wbm is not None:
            if not config['processing']['dem_type'] == 'GETASSE30' and not os.path.isfile(wbm):
                raise FileNotFoundError('External water body mask could not be found: {}'.format(wbm))
        
        dm_path_base = prod_meta["file_base"].format(suffix='dm')
        dm_path = os.path.join(prod_meta['dir_ard_product'], 'annotation', dm_path_base)
        if not os.path.isfile(dm_path):
            log.info(f"creating {os.path.relpath(dm_path, prod_meta['dir_ard_product'])}")
            processor = load_processor(processor_name)
            lsm_encoding = processor.lsm_encoding()
            create_data_mask(outname=dm_path, datasets=sar_assets, extent=extent, epsg=epsg,
                             driver=driver, creation_opt=write_options['dm'],
                             overviews=overviews, overview_resampling=ovr_resampling,
                             dst_nodata=dst_nodata_byte, wbm=wbm,
                             product_type=prod_meta['product_type'],
                             lsm_encoding=lsm_encoding)
        ard_assets['dm'] = dm_path
    
    # create acquisition ID image raster (-id.tif)
    if 'id' in allowed:
        id_path_base = prod_meta["file_base"].format(suffix='id')
        id_path = os.path.join(prod_meta['dir_ard_product'], 'annotation', id_path_base)
        if not os.path.isfile(id_path):
            log.info(f"creating {os.path.relpath(id_path, prod_meta['dir_ard_product'])}")
            create_acq_id_image(outname=id_path, ref_tif=ref_tif,
                                datasets=sar_assets, src_ids=src_ids,
                                extent=extent, epsg=epsg, driver=driver,
                                creation_opt=write_options['id'],
                                overviews=overviews, dst_nodata=dst_nodata_byte)
        ard_assets['id'] = id_path
    
    # create DEM (-em.tif)
    # (if not already converted from processor output)
    if dem_type is not None and 'em' in allowed:
        em_path_base = prod_meta["file_base"].format(suffix='em')
        em_path = os.path.join(prod_meta['dir_ard_product'], 'annotation', em_path_base)
        if not os.path.isfile(em_path):
            log.info(f"creating {os.path.relpath(em_path, prod_meta['dir_ard_product'])}")
            with Raster(ref_tif) as ras:
                tr = ras.res
            log_pyro = logging.getLogger('pyroSAR')
            level = log_pyro.level
            log_pyro.setLevel('NOTSET')
            dem.to_mgrs(dem_type=dem_type, dst=em_path,
                        overviews=overviews, tile=tile, tr=tr,
                        create_options=write_options['em'],
                        pbar=False)
            log_pyro.setLevel(level)
        ard_assets['em'] = em_path
    
    # create color composite VRT (-cc-[gs]-lin.vrt)
    if prod_meta['polarization'] in ['DH', 'DV'] and len(measure_tifs) == 2:
        cc_path = re.sub('[hv]{2}', 'cc', measure_tifs[0]).replace('.tif', '.vrt')
        if not os.path.isfile(cc_path):
            log.info(f"creating {os.path.relpath(cc_path, prod_meta['dir_ard_product'])}")
            create_rgb_vrt(outname=cc_path, infiles=measure_tifs,
                           overviews=overviews, overview_resampling=ovr_resampling)
        key = re.search('cc-[gs]-lin', cc_path).group()
        ard_assets[key] = cc_path
    
    # create log-scaled gamma0|sigma0 nought VRTs (-[vh|vv|hh|hv]-[gs]-log.vrt)
    fun = 'dB'
    args = {'fact': 10}
    scale = None
    for item in measure_tifs:
        target = item.replace('lin.tif', 'log.vrt')
        if not os.path.isfile(target):
            log.info(f"creating {os.path.relpath(target, prod_meta['dir_ard_product'])}")
            create_vrt(src=item, dst=target, fun=fun, scale=scale,
                       args=args, options=vrt_options, overviews=overviews,
                       overview_resampling=ovr_resampling)
        key = re.search('[hv]{2}-[gs]-log', target).group()
        ard_assets[key] = target
    
    # create sigma nought RTC VRTs (-[vh|vv|hh|hv]-s-[lin|log].vrt)
    if 'gs' in allowed:
        gs_path = ard_assets['gs']
        for item in measure_tifs:
            sigma0_rtc_lin = item.replace('g-lin.tif', 's-lin.vrt')
            sigma0_rtc_log = item.replace('g-lin.tif', 's-log.vrt')
            
            if not os.path.isfile(sigma0_rtc_lin):
                log.info(f"creating {os.path.relpath(sigma0_rtc_lin, prod_meta['dir_ard_product'])}")
                create_vrt(src=[item, gs_path], dst=sigma0_rtc_lin, fun='mul',
                           relpaths=True, options=vrt_options, overviews=overviews,
                           overview_resampling=ovr_resampling)
            key = re.search('[hv]{2}-s-lin', sigma0_rtc_lin).group()
            ard_assets[key] = sigma0_rtc_lin
            
            if not os.path.isfile(sigma0_rtc_log):
                log.info(f"creating {os.path.relpath(sigma0_rtc_log, prod_meta['dir_ard_product'])}")
                create_vrt(src=sigma0_rtc_lin, dst=sigma0_rtc_log, fun=fun,
                           scale=scale, options=vrt_options, overviews=overviews,
                           overview_resampling=ovr_resampling, args=args)
            key = key.replace('lin', 'log')
            ard_assets[key] = sigma0_rtc_log
    
    # create gamma nought RTC VRTs (-[vh|vv|hh|hv]-g-[lin|log].vrt)
    if 'sg' in allowed:
        sg_path = ard_assets['sg']
        for item in measure_tifs:
            if not item.endswith('s-lin.tif'):
                continue
            gamma0_rtc_lin = item.replace('s-lin.tif', 'g-lin.vrt')
            gamma0_rtc_log = item.replace('s-lin.tif', 'g-log.vrt')
            
            if not os.path.isfile(gamma0_rtc_lin):
                log.info(f"creating {os.path.relpath(gamma0_rtc_lin, prod_meta['dir_ard_product'])}")
                create_vrt(src=[item, sg_path], dst=gamma0_rtc_lin, fun='mul',
                           relpaths=True, options=vrt_options, overviews=overviews,
                           overview_resampling=ovr_resampling)
            key = re.search('[hv]{2}-g-lin', gamma0_rtc_lin).group()
            ard_assets[key] = gamma0_rtc_lin
            
            if not os.path.isfile(gamma0_rtc_log):
                log.info(f"creating {os.path.relpath(gamma0_rtc_log, prod_meta['dir_ard_product'])}")
                create_vrt(src=gamma0_rtc_lin, dst=gamma0_rtc_log, fun=fun,
                           scale=scale, options=vrt_options, overviews=overviews,
                           overview_resampling=ovr_resampling, args=args)
            key = key.replace('lin', 'log')
            ard_assets[key] = gamma0_rtc_log
    
    # create backscatter wind model (-wm.tif)
    # and wind normalization VRT (-[vv|hh]-s-lin-wn.vrt)
    wm_ref_speed = wm_ref_direction = None
    if 'wm' in annotation:
        wm = []
        wm_ref_speed = []
        wm_ref_direction = []
        for i, ds in enumerate(sar_assets):
            if 'wm' in ds.keys():
                wm.append(ds['wm'])
                # for key in ['wm_ref_speed', 'wm_ref_direction']:
                if 'wm_ref_speed' in ds.keys():
                    wm_ref_speed.append(ds['wm_ref_speed'])
                if 'wm_ref_direction' in ds.keys():
                    wm_ref_direction.append(ds['wm_ref_direction'])
            else:
                scene_base = os.path.basename(src_ids[i].scene)
                raise RuntimeError(f'could not find wind model product for scene {scene_base}')
        
        copol = 'VV' if 'VV' in src_ids[0].polarizations else 'HH'
        copol_sigma0_key = f'{copol.lower()}-s-lin'
        if copol_sigma0_key in ard_assets.keys():
            copol_sigma0 = ard_assets[copol_sigma0_key]
            wn_ard = re.sub(r's-lin\.(?:tif|vrt)', 's-lin-wn.vrt', copol_sigma0)
        else:
            copol_sigma0 = None
            wn_ard = None
        
        outname_base = prod_meta["file_base"].format(suffix='wm')
        wm_ard = os.path.join(prod_meta['dir_ard_product'], 'annotation', outname_base)
        
        gapfill = True if src_ids[0].product == 'GRD' else False
        
        log.info(f"creating {os.path.relpath(wm_ard, prod_meta['dir_ard_product'])}")
        if wn_ard is not None:
            log.info(f"creating {os.path.relpath(wn_ard, prod_meta['dir_ard_product'])}")
        
        wind_normalization(src=wm, dst_wm=wm_ard, dst_wn=wn_ard, measurement=copol_sigma0,
                           gapfill=gapfill, bounds=bounds, epsg=epsg, driver=driver,
                           creation_opt=write_options['wm'],
                           dst_nodata=dst_nodata_float, multithread=multithread)
        ard_assets['wm'] = wm_ard
        ard_assets[f'{copol_sigma0_key}-wn'] = wn_ard
    
    ard_assets = sorted(sorted(list(ard_assets.values()), key=lambda x: os.path.splitext(x)[1]),
                        key=lambda x: os.path.basename(os.path.dirname(x)), reverse=True)
    
    return ard_assets


def append_metadata(
        config: dict[str, dict[str, int | float | str | list[str]]],
        prod_meta: dict[str, str | int | datetime],
        src_ids: list[ID],
        assets: list[str],
        compression: str,
        wm_ref_speed: list[str] | None = None,
        wm_ref_direction: list[str] | None = None
) -> None:
    """
    Append metadata files to an ARD product.
    
    Parameters
    ----------
    config:
        the configuration dictionary
    prod_meta:
        the product metadata as returned by :func:`product_info`
    src_ids:
        the source product objects
    assets:
        a list of assets in the ARD product as returned by :func:`format`.
    compression:
        the used compression algorithm
    wm_ref_speed:
        the wind model reference wind speed files
    wm_ref_direction:
        the wind model reference wind direction files

    Returns
    -------

    """
    meta = meta_dict(config=config, prod_meta=prod_meta,
                     src_ids=src_ids, compression=compression)
    
    append_wind_norm(meta=meta, wm_ref_speed=wm_ref_speed,
                     wm_ref_direction=wm_ref_direction)
    
    # copy support files
    schema_dir = os.path.join(s1ard.__path__[0], 'validation', 'schemas')
    schemas = os.listdir(schema_dir)
    for schema in schemas:
        schema_in = os.path.join(schema_dir, schema)
        schema_out = os.path.join(prod_meta['dir_ard_product'], 'support', schema)
        if not os.path.isfile(schema_out):
            log.info(f"creating {os.path.relpath(schema_out, prod_meta['dir_ard_product'])}")
            shutil.copy(schema_in, schema_out)
    
    if config['metadata']['copy_original']:
        copy_src_meta(ard_dir=prod_meta['dir_ard_product'], src_ids=src_ids)
    if 'OGC' in config['metadata']['format']:
        xml.parse(meta=meta, target=prod_meta['dir_ard_product'],
                  assets=assets, exist_ok=True)
    if 'STAC' in config['metadata']['format']:
        stac.parse(meta=meta, target=prod_meta['dir_ard_product'],
                   assets=assets, exist_ok=True)


def get_datasets(
        scenes: list[str],
        sar_dir: str,
        extent: dict[str, int | float],
        epsg: int,
        processor_name: str
) -> tuple[list[ID], list[dict[str, str]]]:
    """
    Collect processing output for a list of scenes. Reads metadata from all
    source SLC/GRD scenes, finds matching output files in `sar_dir` and
    filters both lists depending on the actual overlap of each SLC/GRD valid
    data coverage with the current MGRS tile geometry. If no output is found
    for any scene the function will raise an error.
    To obtain the extent of valid data coverage, first a binary mask raster
    file is created with the name `datamask.tif`, which is stored in the same
    folder as the processing output as found by :func:`~s1ard.snap.find_datasets`.
    Then, the boundary of this binary mask is computed and stored as `datamask.gpkg`
    (see function :func:`spatialist.vector.boundary`).
    If the provided `extent` does not overlap with this boundary, the output is
    discarded. This scenario might occur when the scene's geometry read from its
    metadata overlaps with the tile but the actual extent of data does not.

    Parameters
    ----------
    scenes:
        List of scenes to process. Either an individual scene or multiple,
        matching scenes (consecutive acquisitions).
    sar_dir:
        The directory containing the SAR datasets processed from the source
        scenes using pyroSAR. The function will raise an error if the processing
        output cannot be found for all scenes in `sar_dir`.
    extent:
        Spatial extent of the MGRS tile, derived from a
        :class:`~spatialist.vector.Vector` object.
    epsg:
        The coordinate reference system as an EPSG code.
    processor_name:
        The name of the used SAR processor. The function `find_datasets` of the
        respective processor module is used.

    Returns
    -------
        List of :class:`~pyroSAR.drivers.ID` objects of all source SLC/GRD scenes
        that overlap with the current MGRS tile and a list of SAR processing output
        files that match each :class:`~pyroSAR.drivers.ID` object of `ids`.
        The format of the latter is a list of dictionaries per scene with keys as
        described by e.g. :func:`s1ard.snap.find_datasets`.
    
    See Also
    --------
    :func:`s1ard.snap.find_datasets`
    """
    processor = load_processor(processor_name)
    ids = identify_many(scenes, sortkey='start')
    datasets = []
    for i, _id in enumerate(ids):
        log.debug(f'collecting processing output for scene {os.path.basename(_id.scene)}')
        files = processor.find_datasets(scene=_id.scene, outdir=sar_dir, epsg=epsg)
        if files is not None:
            base = os.path.splitext(os.path.basename(_id.scene))[0]
            ocn = re.sub('(?:SLC_|GRD[FHM])_1', 'OCN__2', base)[:-5]
            # allow 1 second tolerance
            s_start = int(ocn[31])
            s_stop = int(ocn[47])
            ocn_list = list(ocn)
            s = 1
            ocn_list[31] = f'[{s_start - s}{s_start}{s_start + s}]'
            ocn_list[47] = f'[{s_stop - s}{s_stop}{s_stop + s}]'
            ocn = ''.join(ocn_list)
            log.debug(f'searching for OCN products with pattern {ocn}')
            ocn_match = finder(target=sar_dir, matchlist=[ocn], regex=True,
                               foldermode=2, recursive=False)
            if len(ocn_match) > 0:
                for v in ['owiNrcsCmod', 'owiEcmwfWindSpeed', 'owiEcmwfWindDirection']:
                    ocn_tif = os.path.join(ocn_match[0], f'{v}.tif')
                    if os.path.isfile(ocn_tif):
                        if v.endswith('Speed'):
                            files['wm_ref_speed'] = ocn_tif
                        elif v.endswith('Direction'):
                            files['wm_ref_direction'] = ocn_tif
                        else:
                            files['wm'] = ocn_tif
            datasets.append(files)
        else:
            base = os.path.basename(_id.scene)
            msg = f'cannot find processing output for scene {base} and CRS EPSG:{epsg}'
            raise RuntimeError(msg)
    
    i = 0
    while i < len(datasets):
        log.debug(f'checking tile overlap for scene {os.path.basename(ids[i].scene)}')
        measurements = [datasets[i][x] for x in datasets[i].keys()
                        if re.search('[gs]-lin', x)]
        dm_ras = os.path.join(os.path.dirname(measurements[0]), 'datamask.tif')
        dm_vec = dm_ras.replace('.tif', '.gpkg')
        dm_vec = datamask(measurement=measurements[0], dm_ras=dm_ras, dm_vec=dm_vec)
        if dm_vec is None:
            del ids[i], datasets[i]
            continue
        with Lock(dm_vec, soft=True):
            with Vector(dm_vec) as bounds:
                with bbox(extent, epsg) as tile_geom:
                    inter = intersect(bounds, tile_geom)
                    if inter is not None:
                        with Raster(dm_ras) as ras:
                            inter_min = ras.res[0] * ras.res[1]
                        if inter.getArea() < inter_min:
                            inter.close()
                            inter = None
                    if inter is None:
                        log.debug('no overlap, removing scene')
                        del ids[i]
                        del datasets[i]
                    else:
                        log.debug('overlap detected')
                        # Add dm_ras to the datasets if it overlaps with the current tile
                        datasets[i]['datamask'] = dm_ras
                        i += 1
                        inter.close()
    return ids, datasets


def wind_normalization(
        src: list[str],
        dst_wm: str,
        dst_wn: str | None,
        measurement: str | None,
        gapfill: bool,
        bounds: list[float],
        epsg: int,
        driver: str,
        creation_opt: list[str],
        dst_nodata: float,
        multithread: bool,
        resolution: int = 915
) -> None:
    """
    Create ORB wind normalization layers. A wind model annotation layer is
    created and optionally a wind normalization VRT.
    
    Parameters
    ----------
    src:
        A list of OCN products as prepared by :func:`s1ard.ocn.extract`
    dst_wm:
        The name of the wind model layer in the ARD product
    dst_wn:
        The name of the wind normalization VRT. If None, no VRT will be created.
        Requires `measurement` to point to a file.
    measurement:
        The name of the measurement file used for wind normalization in `dst_wn`.
        If None, no wind normalization VRT will be created.
    gapfill:
        Perform additional gap filling (:func:`s1ard.ocn.gapfill`)?
        This is recommended if the Level-1 source product of `measurement` is GRD
        in which case gaps are introduced between subsequently acquired scenes.
    bounds:
        the bounds of the MGRS tile
    epsg:
        The EPSG code of the MGRS tile
    driver:
        GDAL driver to use for raster file creation.
    creation_opt:
        GDAL creation options to use for raster file creation. Should match
        specified GDAL driver.
    dst_nodata:
        Nodata value to write to the output raster.
    multithread:
        Should `gdalwarp` use multithreading?
    resolution:
        The target pixel resolution in meters. 915 is chosen as default because
        it is closest to the OCN product resolution (1000) and still fits into
        the MGRS bounds (``109800 % 915 == 0``).
    
    Returns
    -------
    
    """
    if len(src) > 1:
        cmod_mosaic = get_tmp_name(suffix='.tif')
        gdalwarp(src=src, dst=cmod_mosaic)
        if gapfill:
            cmod_geo = get_tmp_name(suffix='.tif')
            ocn.gapfill(src=cmod_mosaic, dst=cmod_geo, md=2, si=1)
            os.remove(cmod_mosaic)
        else:
            cmod_geo = cmod_mosaic
        cmod_geo_tmp = True
    else:
        cmod_geo = src[0]
        cmod_geo_tmp = False
    
    if not os.path.isfile(dst_wm):
        log.info(f'creating {dst_wm}')
        gdalwarp(src=cmod_geo,
                 dst=dst_wm,
                 outputBounds=bounds,
                 dstSRS=f'EPSG:{epsg}',
                 xRes=resolution, yRes=resolution,
                 resampleAlg='bilinear',
                 format=driver,
                 dstNodata=dst_nodata,
                 multithread=multithread,
                 creationOptions=creation_opt)
    
    if cmod_geo_tmp:
        os.remove(cmod_geo)
    
    if dst_wn is not None and measurement is not None:
        if not os.path.isfile(dst_wn):
            log.info(f'creating {dst_wn}')
            with Raster(measurement) as ras:
                xres, yres = ras.res
            create_vrt(src=[measurement, dst_wm], dst=dst_wn, fun='div',
                       options={'xRes': xres, 'yRes': yres}, relpaths=True)
