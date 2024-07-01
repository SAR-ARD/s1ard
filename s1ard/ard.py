import os
import re
import time
import shutil
import tempfile
from datetime import datetime, timezone
import numpy as np
from lxml import etree
from time import gmtime, strftime
from copy import deepcopy
from scipy.interpolate import griddata
from osgeo import gdal
from spatialist.vector import Vector, vectorize, boundary, bbox, intersect
from spatialist.raster import Raster, rasterize, Dtype
from spatialist.auxil import gdalwarp, gdalbuildvrt
from spatialist.ancillary import finder
from pyroSAR import identify, identify_many
import s1ard
from s1ard import dem, ocn
from s1ard.metadata import extract, xml, stac
from s1ard.metadata.mapping import LERC_ERR_THRES
from s1ard.ancillary import generate_unique_id, vrt_add_overviews
from s1ard.metadata.extract import copy_src_meta, get_src_meta, find_in_annotation
from s1ard.snap import find_datasets
import logging

log = logging.getLogger('s1ard')


def format(config, product_type, scenes, datadir, outdir, tile, extent, epsg, wbm=None,
           dem_type=None, multithread=True, compress=None, overviews=None,
           annotation=None, update=False):
    """
    Finalizes the generation of Sentinel-1 Analysis Ready Data (ARD) products after SAR processing has finished.
    This includes the following:
    
    - Creating all measurement and annotation datasets in Cloud Optimized GeoTIFF (COG) format
    - Creating additional annotation datasets in Virtual Raster Tile (VRT) format
    - Applying the ARD product directory structure & naming convention
    - Generating metadata in XML and JSON formats for the ARD product as well as source SLC datasets
    
    Parameters
    ----------
    config: dict
        Dictionary of the parsed config parameters for the current process.
    product_type: str
        The type of ARD product to be generated. Options: 'NRB' or 'ORB'.
    scenes: list[str]
        List of scenes to process. Either a single scene or multiple, matching scenes (consecutive acquisitions).
        All scenes are expected to overlap with `extent` and an error will be thrown if the processing output
        cannot be found for any of the scenes.
    datadir: str
        The directory containing the SAR datasets processed from the source scenes using pyroSAR.
    outdir: str
        The directory to write the final files to.
    tile: str
        ID of an MGRS tile.
    extent: dict
        Spatial extent of the MGRS tile, derived from a :class:`~spatialist.vector.Vector` object.
    epsg: int
        The CRS used for the ARD product; provided as an EPSG code.
    wbm: str or None
        Path to a water body mask file with the dimensions of an MGRS tile.
    dem_type: str or None
        if defined, a DEM layer will be added to the product. The suffix `em` (elevation model) is used.
        Default `None`: do not add a DEM layer.
    multithread: bool
        Should `gdalwarp` use multithreading? Default is True. The number of threads used, can be adjusted in the
        `config.ini` file with the parameter `gdal_threads`.
    compress: str or None
        Compression algorithm to use. See https://gdal.org/drivers/raster/gtiff.html#creation-options for options.
        Defaults to 'LERC_DEFLATE'.
    overviews: list[int] or None
        Internal overview levels to be created for each GeoTIFF file. Defaults to [2, 4, 9, 18, 36]
    annotation: list[str] or None
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
    update: bool
        modify existing products so that only missing files are re-created?
    
    Returns
    -------
    str
        Either the time spent executing the function in seconds or 'Already processed - Skip!'
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
    
    # determine processing timestamp and generate unique ID
    start_time = time.time()
    proc_time = datetime.now(timezone.utc)
    t = proc_time.isoformat().encode()
    product_id = generate_unique_id(encoded_str=t)
    
    src_ids, datasets_sar = get_datasets(scenes=scenes, datadir=datadir, extent=extent, epsg=epsg)
    if len(src_ids) == 0:
        log.error(f'None of the processed scenes overlap with the current tile {tile}')
        return
    
    if annotation is not None:
        allowed = []
        for key in datasets_sar[0]:
            c1 = re.search('[gs]-lin', key)
            c2 = key in annotation
            c3 = key.startswith('np') and 'np' in annotation
            if c1 or c2 or c3:
                allowed.append(key)
    else:
        allowed = [key for key in datasets_sar[0].keys() if re.search('[gs]-lin', key)]
        annotation = []
    for item in ['em', 'id']:
        if item in annotation:
            allowed.append(item)
    
    # GDAL output bounds
    bounds = [extent['xmin'], extent['ymin'], extent['xmax'], extent['ymax']]
    
    ard_start, ard_stop = calc_product_start_stop(src_ids=src_ids, extent=extent, epsg=epsg)
    meta = {'mission': src_ids[0].sensor,
            'mode': src_ids[0].meta['acquisition_mode'],
            'ard_spec': product_type,
            'polarization': {"['HH']": 'SH',
                             "['VV']": 'SV',
                             "['HH', 'HV']": 'DH',
                             "['VV', 'VH']": 'DV'}[str(src_ids[0].polarizations)],
            'start': ard_start,
            'orbitnumber': src_ids[0].meta['orbitNumbers_abs']['start'],
            'datatake': hex(src_ids[0].meta['frameNumber']).replace('x', '').upper(),
            'tile': tile,
            'id': product_id}
    meta_lower = dict((k, v.lower() if not isinstance(v, int) else v) for k, v in meta.items())
    skeleton_dir = '{mission}_{mode}_{ard_spec}__1S{polarization}_{start}_{orbitnumber:06}_{datatake:0>6}_{tile}_{id}'
    skeleton_files = '{mission}-{mode}-{ard_spec}-{start}-{orbitnumber:06}-{datatake:0>6}-{tile}-{suffix}.tif'
    
    ard_base = skeleton_dir.format(**meta)
    log.info(f'product name: {os.path.join(outdir, ard_base)}')
    existing = finder(outdir, [ard_base.replace(product_id, '*')], foldermode=2)
    if len(existing) > 0:
        if not update:
            return 'Already processed - Skip!'
        else:
            ard_dir = existing[0]
    else:
        ard_dir = os.path.join(outdir, ard_base)
    os.makedirs(ard_dir, exist_ok=True)
    subdirectories = ['measurement', 'annotation', 'source', 'support']
    for subdirectory in subdirectories:
        os.makedirs(os.path.join(ard_dir, subdirectory), exist_ok=True)
    
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
    datasets_ard = dict()
    for key in list(datasets_sar[0].keys()):
        if key in ['dm', 'wm'] or key not in LERC_ERR_THRES.keys() or key not in allowed:
            # raster files for keys 'dm' and 'wm' are created later
            continue
        
        meta_lower['suffix'] = key
        outname_base = skeleton_files.format(**meta_lower)
        if re.search('[gs]-lin', key):
            subdir = 'measurement'
        else:
            subdir = 'annotation'
        outname = os.path.join(ard_dir, subdir, outname_base)
        
        if not os.path.isfile(outname):
            log.info(f'creating {outname}')
            images = [ds[key] for ds in datasets_sar]
            ras = None
            if len(images) > 1:
                ras = Raster(images, list_separate=False)
                source = ras.filename
            else:
                source = tempfile.NamedTemporaryFile(suffix='.vrt').name
                gdalbuildvrt(src=images[0], dst=source)
            
            # modify temporary VRT to make sure overview levels and resampling are properly applied
            vrt_add_overviews(vrt=source, overviews=overviews, resampling=ovr_resampling)
            
            options = {'format': driver, 'outputBounds': bounds,
                       'dstNodata': dst_nodata_float, 'multithread': multithread,
                       'creationOptions': write_options[key]}
            
            gdalwarp(src=source, dst=outname, **options)
            if ras is not None:
                ras.close()
        datasets_ard[key] = outname
    
    # define a reference raster from the annotation datasets and list all gamma0/sigma0 backscatter measurement rasters
    measure_tifs = [v for k, v in datasets_ard.items() if re.search('[gs]-lin', k)]
    ref_key = list(datasets_ard.keys())[0]
    ref_tif = datasets_ard[ref_key]
    
    # create data mask raster (-dm.tif)
    if 'dm' in allowed:
        if wbm is not None:
            if not config['dem_type'] == 'GETASSE30' and not os.path.isfile(wbm):
                raise FileNotFoundError('External water body mask could not be found: {}'.format(wbm))
        
        dm_path = ref_tif.replace(f'-{ref_key}.tif', '-dm.tif')
        if not os.path.isfile(dm_path):
            log.info(f'creating {dm_path}')
            create_data_mask(outname=dm_path, datasets=datasets_sar, extent=extent, epsg=epsg,
                             driver=driver, creation_opt=write_options['dm'],
                             overviews=overviews, overview_resampling=ovr_resampling,
                             dst_nodata=dst_nodata_byte, wbm=wbm, product_type=product_type)
        datasets_ard['dm'] = dm_path
    
    # create acquisition ID image raster (-id.tif)
    if 'id' in allowed:
        id_path = ref_tif.replace(f'-{ref_key}.tif', '-id.tif')
        if not os.path.isfile(id_path):
            log.info(f'creating {id_path}')
            create_acq_id_image(outname=id_path, ref_tif=ref_tif,
                                datasets=datasets_sar, src_ids=src_ids,
                                extent=extent, epsg=epsg, driver=driver,
                                creation_opt=write_options['id'],
                                overviews=overviews, dst_nodata=dst_nodata_byte)
        datasets_ard['id'] = id_path
    
    # create DEM (-em.tif)
    if dem_type is not None and 'em' in allowed:
        em_path = ref_tif.replace(f'-{ref_key}.tif', '-em.tif')
        if not os.path.isfile(em_path):
            log.info(f'creating {em_path}')
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
        datasets_ard['em'] = em_path
    
    # create color composite VRT (-cc-[gs]-lin.vrt)
    if meta['polarization'] in ['DH', 'DV'] and len(measure_tifs) == 2:
        cc_path = re.sub('[hv]{2}', 'cc', measure_tifs[0]).replace('.tif', '.vrt')
        if not os.path.isfile(cc_path):
            log.info(f'creating {cc_path}')
            create_rgb_vrt(outname=cc_path, infiles=measure_tifs,
                           overviews=overviews, overview_resampling=ovr_resampling)
        key = re.search('cc-[gs]-lin', cc_path).group()
        datasets_ard[key] = cc_path
    
    # create log-scaled gamma0|sigma0 nought VRTs (-[vh|vv|hh|hv]-[gs]-log.vrt)
    fun = 'dB'
    args = {'fact': 10}
    scale = None
    for item in measure_tifs:
        target = item.replace('lin.tif', 'log.vrt')
        if not os.path.isfile(target):
            log.info(f'creating {target}')
            create_vrt(src=item, dst=target, fun=fun, scale=scale,
                       args=args, options=vrt_options, overviews=overviews,
                       overview_resampling=ovr_resampling)
        key = re.search('[hv]{2}-[gs]-log', target).group()
        datasets_ard[key] = target
    
    # create sigma nought RTC VRTs (-[vh|vv|hh|hv]-s-[lin|log].vrt)
    if 'gs' in allowed:
        gs_path = datasets_ard['gs']
        for item in measure_tifs:
            sigma0_rtc_lin = item.replace('g-lin.tif', 's-lin.vrt')
            sigma0_rtc_log = item.replace('g-lin.tif', 's-log.vrt')
            
            if not os.path.isfile(sigma0_rtc_lin):
                log.info(f'creating {sigma0_rtc_lin}')
                create_vrt(src=[item, gs_path], dst=sigma0_rtc_lin, fun='mul',
                           relpaths=True, options=vrt_options, overviews=overviews,
                           overview_resampling=ovr_resampling)
            key = re.search('[hv]{2}-s-lin', sigma0_rtc_lin).group()
            datasets_ard[key] = sigma0_rtc_lin
            
            if not os.path.isfile(sigma0_rtc_log):
                log.info(f'creating {sigma0_rtc_log}')
                create_vrt(src=sigma0_rtc_lin, dst=sigma0_rtc_log, fun=fun,
                           scale=scale, options=vrt_options, overviews=overviews,
                           overview_resampling=ovr_resampling, args=args)
            key = key.replace('lin', 'log')
            datasets_ard[key] = sigma0_rtc_log
    
    # create gamma nought RTC VRTs (-[vh|vv|hh|hv]-g-[lin|log].vrt)
    if 'sg' in allowed:
        sg_path = datasets_ard['sg']
        for item in measure_tifs:
            if not item.endswith('s-lin.tif'):
                continue
            gamma0_rtc_lin = item.replace('s-lin.tif', 'g-lin.vrt')
            gamma0_rtc_log = item.replace('s-lin.tif', 'g-log.vrt')
            
            if not os.path.isfile(gamma0_rtc_lin):
                log.info(f'creating {gamma0_rtc_lin}')
                create_vrt(src=[item, sg_path], dst=gamma0_rtc_lin, fun='mul',
                           relpaths=True, options=vrt_options, overviews=overviews,
                           overview_resampling=ovr_resampling)
            key = re.search('[hv]{2}-g-lin', gamma0_rtc_lin).group()
            datasets_ard[key] = gamma0_rtc_lin
            
            if not os.path.isfile(gamma0_rtc_log):
                log.info(f'creating {gamma0_rtc_log}')
                create_vrt(src=gamma0_rtc_lin, dst=gamma0_rtc_log, fun=fun,
                           scale=scale, options=vrt_options, overviews=overviews,
                           overview_resampling=ovr_resampling, args=args)
            key = key.replace('lin', 'log')
            datasets_ard[key] = gamma0_rtc_log
    
    # create backscatter wind model (-wm.tif)
    # and wind normalization VRT (-[vv|hh]-s-lin-wn.vrt)
    wm_ref_files = None
    if 'wm' in annotation:
        wm = []
        wm_ref_files = []
        for i, ds in enumerate(datasets_sar):
            if 'wm' in ds.keys():
                wm.append(ds['wm'])
                for key in ['wm_ref_speed', 'wm_ref_direction']:
                    if key in ds.keys():
                        wm_ref_files.append(ds[key])
            else:
                scene_base = os.path.basename(src_ids[i].scene)
                raise RuntimeError(f'could not find wind model product for scene {scene_base}')
        
        copol = 'VV' if 'VV' in src_ids[0].polarizations else 'HH'
        copol_sigma0_key = f'{copol.lower()}-s-lin'
        if copol_sigma0_key in datasets_ard.keys():
            copol_sigma0 = datasets_ard[copol_sigma0_key]
            wn_ard = re.sub(r's-lin\.(?:tif|vrt)', 's-lin-wn.vrt', copol_sigma0)
        else:
            copol_sigma0 = None
            wn_ard = None
        
        meta_lower['suffix'] = 'wm'
        outname_base = skeleton_files.format(**meta_lower)
        wm_ard = os.path.join(ard_dir, 'annotation', outname_base)
        
        gapfill = True if src_ids[0].product == 'GRD' else False
        
        wind_normalization(src=wm, dst_wm=wm_ard, dst_wn=wn_ard, measurement=copol_sigma0,
                           gapfill=gapfill, bounds=bounds, epsg=epsg, driver=driver,
                           creation_opt=write_options['wm'],
                           dst_nodata=dst_nodata_float, multithread=multithread)
        datasets_ard['wm'] = wm_ard
        datasets_ard[f'{copol_sigma0_key}-wn'] = wn_ard
    
    # copy support files
    schema_dir = os.path.join(s1ard.__path__[0], 'validation', 'schemas')
    schemas = os.listdir(schema_dir)
    for schema in schemas:
        schema_in = os.path.join(schema_dir, schema)
        schema_out = os.path.join(ard_dir, 'support', schema)
        if not os.path.isfile(schema_out):
            log.info(f'creating {schema_out}')
            shutil.copy(schema_in, schema_out)
    
    # create metadata files in XML and (STAC) JSON formats
    start = datetime.strptime(ard_start, '%Y%m%dT%H%M%S')
    stop = datetime.strptime(ard_stop, '%Y%m%dT%H%M%S')
    meta = extract.meta_dict(config=config, target=ard_dir, src_ids=src_ids, sar_dir=datadir,
                             proc_time=proc_time, start=start, stop=stop, compression=compress,
                             product_type=product_type, wm_ref_files=wm_ref_files)
    ard_assets = sorted(sorted(list(datasets_ard.values()), key=lambda x: os.path.splitext(x)[1]),
                        key=lambda x: os.path.basename(os.path.dirname(x)), reverse=True)
    if config['meta']['copy_original']:
        copy_src_meta(ard_dir=ard_dir, src_ids=src_ids)
    if 'OGC' in config['meta']['format']:
        xml.parse(meta=meta, target=ard_dir, assets=ard_assets, exist_ok=True)
    if 'STAC' in config['meta']['format']:
        stac.parse(meta=meta, target=ard_dir, assets=ard_assets, exist_ok=True)
    return str(round((time.time() - start_time), 2))


def get_datasets(scenes, datadir, extent, epsg):
    """
    Collect processing output for a list of scenes.
    Reads metadata from all source SLC/GRD scenes, finds matching output files in `datadir`
    and filters both lists depending on the actual overlap of each SLC/GRD valid data coverage
    with the current MGRS tile geometry. If no output is found for any scene the function will raise an error.
    To obtain the extent of valid data coverage, first a binary
    mask raster file is created with the name `datamask.tif`, which is stored in the same folder as
    the processing output as found by :func:`~s1ard.snap.find_datasets`. Then, the boundary of this
    binary mask is computed and stored as `datamask.gpkg` (see function :func:`spatialist.vector.boundary`).
    If the provided `extent` does not overlap with this boundary, the output is discarded. This scenario
    might occur when the scene's geometry read from its metadata overlaps with the tile but the actual
    extent of data does not.

    Parameters
    ----------
    scenes: list[str]
        List of scenes to process. Either an individual scene or multiple, matching scenes (consecutive acquisitions).
    datadir: str
        The directory containing the SAR datasets processed from the source scenes using pyroSAR.
        The function will raise an error if the processing output cannot be found for all scenes in `datadir`.
    extent: dict
        Spatial extent of the MGRS tile, derived from a :class:`~spatialist.vector.Vector` object.
    epsg: int
        The coordinate reference system as an EPSG code.

    Returns
    -------
    ids: list[:class:`pyroSAR.drivers.ID`]
        List of :class:`~pyroSAR.drivers.ID` objects of all source SLC/GRD scenes that overlap with the current MGRS tile.
    datasets: list[dict]
        List of SAR processing output files that match each :class:`~pyroSAR.drivers.ID` object of `ids`.
        The format is a list of dictionaries per scene with keys as described by e.g. :func:`s1ard.snap.find_datasets`.
    
    See Also
    --------
    :func:`s1ard.snap.find_datasets`
    """
    ids = identify_many(scenes)
    datasets = []
    for i, _id in enumerate(ids):
        files = find_datasets(scene=_id.scene, outdir=datadir, epsg=epsg)
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
            ocn_match = finder(target=datadir, matchlist=[ocn], regex=True, foldermode=2)
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
            raise RuntimeError(f'cannot find processing output for scene {base} and CRS EPSG:{epsg}')
    
    i = 0
    while i < len(datasets):
        measurements = [datasets[i][x] for x in datasets[i].keys() if re.search('[gs]-lin', x)]
        dm_ras = os.path.join(os.path.dirname(measurements[0]), 'datamask.tif')
        dm_vec = dm_ras.replace('.tif', '.gpkg')
        
        if not os.path.isfile(dm_ras):
            with Raster(measurements[0]) as ras:
                arr = ras.array()
                mask = ~np.isnan(arr)
                del arr
                # remove scene if file does not contain valid data
                if len(mask[mask == 1]) == 0:
                    del ids[i], datasets[i]
                    continue
                with vectorize(target=mask, reference=ras) as vec:
                    with boundary(vec, expression="value=1") as bounds:
                        if not os.path.isfile(dm_ras):
                            rasterize(vectorobject=bounds, reference=ras, outname=dm_ras)
                        if not os.path.isfile(dm_vec):
                            bounds.write(outfile=dm_vec)
                del mask
        if not os.path.isfile(dm_vec):
            with Raster(dm_ras) as ras:
                mask = ras.array().astype('bool')
                # remove scene if file does not contain valid data
                if len(mask[mask == 1]) == 0:
                    del ids[i], datasets[i]
                    continue
                with vectorize(target=mask, reference=ras) as vec:
                    boundary(vec, expression="value=1", outname=dm_vec)
                del mask
        with Vector(dm_vec) as bounds:
            with bbox(extent, epsg) as tile_geom:
                inter = intersect(bounds, tile_geom)
                if inter is None:
                    del ids[i]
                    del datasets[i]
                else:
                    # Add dm_ras to the datasets if it overlaps with the current tile
                    datasets[i]['datamask'] = dm_ras
                    i += 1
                    inter.close()
    return ids, datasets


def create_vrt(src, dst, fun, relpaths=False, scale=None, offset=None, dtype=None,
               args=None, options=None, overviews=None, overview_resampling=None):
    """
    Creates a VRT file for the specified source dataset(s) and adds a pixel function that should be applied on the fly
    when opening the VRT file.

    Parameters
    ----------
    src: str or list[str]
        The input dataset(s).
    dst: str
        The output dataset.
    fun: str
        A `PixelFunctionType` that should be applied on the fly when opening the VRT file. The function is applied to a
        band that derives its pixel information from the source bands. A list of possible options can be found here:
        https://gdal.org/drivers/raster/vrt.html#default-pixel-functions.
        Furthermore, the option 'decibel' can be specified, which will implement a custom pixel function that uses
        Python code for decibel conversion (10*log10).
    relpaths: bool
        Should all `SourceFilename` XML elements with attribute `@relativeToVRT="0"` be updated to be paths relative to
        the output VRT file? Default is False.
    scale: int or None
         The scale that should be applied when computing “real” pixel values from scaled pixel values on a raster band.
         Will be ignored if `fun='decibel'`.
    offset: float or None
        The offset that should be applied when computing “real” pixel values from scaled pixel values on a raster band.
        Will be ignored if `fun='decibel'`.
    dtype: str or None
        the data type of the written VRT file; default None: same data type as source data.
        data type notations of GDAL (e.g. `Float32`) and numpy (e.g. `int8`) are supported.
    args: dict or None
        arguments for `fun` passed as `PixelFunctionArguments`. Requires GDAL>=3.5 to be read.
    options: dict or None
        Additional parameters passed to `gdal.BuildVRT`.
    overviews: list[int] or None
        Internal overview levels to be created for each raster file.
    overview_resampling: str or None
        Resampling method for overview levels.

    Examples
    --------
    linear gamma0 backscatter as input:

    >>> src = 's1a-iw-nrb-20220601t052704-043465-0530a1-32tpt-vh-g-lin.tif'

    decibel scaling I:
    use `log10` pixel function and additional `Scale` parameter.
    Known to display well in QGIS, but `Scale` is ignored when reading array in Python.

    >>> dst = src.replace('-lin.tif', '-log1.vrt')
    >>> create_vrt(src=src, dst=dst, fun='log10', scale=10)

    decibel scaling II:
    use custom Python pixel function. Requires additional environment variable GDAL_VRT_ENABLE_PYTHON set to YES.

    >>> dst = src.replace('-lin.tif', '-log2.vrt')
    >>> create_vrt(src=src, dst=dst, fun='decibel')

    decibel scaling III:
    use `dB` pixel function with additional `PixelFunctionArguments`. Works best but requires GDAL>=3.5.

    >>> dst = src.replace('-lin.tif', '-log3.vrt')
    >>> create_vrt(src=src, dst=dst, fun='dB', args={'fact': 10})
    """
    options = {} if options is None else options
    gdalbuildvrt(src=src, dst=dst, **options)
    tree = etree.parse(dst)
    root = tree.getroot()
    band = tree.find('VRTRasterBand')
    band.attrib['subClass'] = 'VRTDerivedRasterBand'
    
    if dtype is not None:
        band.attrib['dataType'] = Dtype(dtype).gdalstr
    
    if fun == 'decibel':
        pxfun_language = etree.SubElement(band, 'PixelFunctionLanguage')
        pxfun_language.text = 'Python'
        pxfun_type = etree.SubElement(band, 'PixelFunctionType')
        pxfun_type.text = fun
        pxfun_code = etree.SubElement(band, 'PixelFunctionCode')
        pxfun_code.text = etree.CDATA("""
    import numpy as np
    def decibel(in_ar, out_ar, xoff, yoff, xsize, ysize, raster_xsize, raster_ysize, buf_radius, gt, **kwargs):
        np.multiply(np.log10(in_ar[0], where=in_ar[0]>0.0, out=out_ar, dtype='float32'), 10.0, out=out_ar, dtype='float32')
        """)
    else:
        pixfun_type = etree.SubElement(band, 'PixelFunctionType')
        pixfun_type.text = fun
        if args is not None:
            arg = etree.SubElement(band, 'PixelFunctionArguments')
            for key, value in args.items():
                arg.attrib[key] = str(value)
        if scale is not None:
            sc = etree.SubElement(band, 'Scale')
            sc.text = str(scale)
        if offset is not None:
            off = etree.SubElement(band, 'Offset')
            off.text = str(offset)
    
    if any([overviews, overview_resampling]) is not None:
        ovr = tree.find('OverviewList')
        if ovr is None:
            ovr = etree.SubElement(root, 'OverviewList')
        if overview_resampling is not None:
            ovr.attrib['resampling'] = overview_resampling.lower()
        if overviews is not None:
            ov = str(overviews)
            for x in ['[', ']', ',']:
                ov = ov.replace(x, '')
            ovr.text = ov
    
    if relpaths:
        srcfiles = tree.xpath('//SourceFilename[@relativeToVRT="0"]')
        for srcfile in srcfiles:
            repl = os.path.relpath(srcfile.text, start=os.path.dirname(dst))
            repl = repl.replace('\\', '/')
            srcfile.text = repl
            srcfile.attrib['relativeToVRT'] = '1'
    
    etree.indent(root)
    tree.write(dst, pretty_print=True, xml_declaration=False, encoding='utf-8')


def create_rgb_vrt(outname, infiles, overviews, overview_resampling):
    """
    Creation of the color composite VRT file.

    Parameters
    ----------
    outname: str
        Full path to the output VRT file.
    infiles: list[str]
        A list of paths pointing to the linear scaled measurement backscatter files.
    overviews: list[int]
        Internal overview levels to be defined for the created VRT file.
    overview_resampling: str
        Resampling method applied to overview pyramids.
    """
    
    # make sure order is right and co-polarization (VV or HH) is first
    pols = [re.search('[hv]{2}', os.path.basename(f)).group() for f in infiles]
    if pols[1] in ['vv', 'hh']:
        infiles.reverse()
        pols.reverse()
    
    # format overview levels
    ov = str(overviews)
    for x in ['[', ']', ',']:
        ov = ov.replace(x, '')
    
    # create VRT file and change its content
    gdalbuildvrt(src=infiles, dst=outname, separate=True)
    
    tree = etree.parse(outname)
    root = tree.getroot()
    srs = tree.find('SRS').text
    geotrans = tree.find('GeoTransform').text
    bands = tree.findall('VRTRasterBand')
    vrt_nodata = bands[0].find('NoDataValue').text
    complex_src = [band.find('ComplexSource') for band in bands]
    for cs in complex_src:
        cs.remove(cs.find('NODATA'))
    
    new_band = etree.SubElement(root, 'VRTRasterBand',
                                attrib={'dataType': 'Float32', 'band': '3',
                                        'subClass': 'VRTDerivedRasterBand'})
    new_band_na = etree.SubElement(new_band, 'NoDataValue')
    new_band_na.text = 'nan'
    pxfun_type = etree.SubElement(new_band, 'PixelFunctionType')
    pxfun_type.text = 'mul'
    for cs in complex_src:
        new_band.append(deepcopy(cs))
    
    src = new_band.findall('ComplexSource')[1]
    fname = src.find('SourceFilename')
    fname_old = fname.text
    src_attr = src.find('SourceProperties').attrib
    fname.text = etree.CDATA("""
    <VRTDataset rasterXSize="{rasterxsize}" rasterYSize="{rasterysize}">
        <SRS dataAxisToSRSAxisMapping="1,2">{srs}</SRS>
        <GeoTransform>{geotrans}</GeoTransform>
        <VRTRasterBand dataType="{dtype}" band="1" subClass="VRTDerivedRasterBand">
            <NoDataValue>{vrt_nodata}</NoDataValue>
            <PixelFunctionType>{px_fun}</PixelFunctionType>
            <ComplexSource>
              <SourceFilename relativeToVRT="1">{fname}</SourceFilename>
              <SourceBand>1</SourceBand>
              <SourceProperties RasterXSize="{rasterxsize}" RasterYSize="{rasterysize}" DataType="{dtype}" BlockXSize="{blockxsize}" BlockYSize="{blockysize}"/>
              <SrcRect xOff="0" yOff="0" xSize="{rasterxsize}" ySize="{rasterysize}"/>
              <DstRect xOff="0" yOff="0" xSize="{rasterxsize}" ySize="{rasterysize}"/>
            </ComplexSource>
        </VRTRasterBand>
        <OverviewList resampling="{ov_resampling}">{ov}</OverviewList>
    </VRTDataset>
    """.format(rasterxsize=src_attr['RasterXSize'], rasterysize=src_attr['RasterYSize'], srs=srs, geotrans=geotrans,
               dtype=src_attr['DataType'], px_fun='inv', fname=fname_old, vrt_nodata=vrt_nodata,
               blockxsize=src_attr['BlockXSize'], blockysize=src_attr['BlockYSize'],
               ov_resampling=overview_resampling.lower(), ov=ov))
    
    bands = tree.findall('VRTRasterBand')
    for band, col in zip(bands, ['Red', 'Green', 'Blue']):
        color = etree.Element('ColorInterp')
        color.text = col
        band.insert(0, color)
    
    ovr = etree.SubElement(root, 'OverviewList', attrib={'resampling': overview_resampling.lower()})
    ovr.text = ov
    
    etree.indent(root)
    tree.write(outname, pretty_print=True, xml_declaration=False, encoding='utf-8')


def calc_product_start_stop(src_ids, extent, epsg):
    """
    Calculates the start and stop times of the ARD product.
    The geolocation grid points including their azimuth time information are extracted first from the metadata of each
    source SLC. These grid points are then used to interpolate the azimuth time for the lower right and upper left
    (ascending) or upper right and lower left (descending) corners of the MGRS tile extent.

    Parameters
    ----------
    src_ids: list[pyroSAR.drivers.ID]
        List of :class:`~pyroSAR.drivers.ID` objects of all source SLC scenes that overlap with the current MGRS tile.
    extent: dict
        Spatial extent of the MGRS tile, derived from a :class:`~spatialist.vector.Vector` object.
    epsg: int
        The coordinate reference system as an EPSG code.

    Returns
    -------
    tuple[str]
        Start and stop time of the ARD product formatted as `YYYYmmddTHHMMSS` in UTC.
    """
    with bbox(extent, epsg) as tile_geom:
        tile_geom.reproject(4326)
        ext = tile_geom.extent
        ul = (ext['xmin'], ext['ymax'])
        ur = (ext['xmax'], ext['ymax'])
        lr = (ext['xmax'], ext['ymin'])
        ll = (ext['xmin'], ext['ymin'])
        tile_geom = None
    
    slc_dict = {}
    for i, sid in enumerate(src_ids):
        uid = os.path.basename(sid.scene).split('.')[0][-4:]
        slc_dict[uid] = get_src_meta(sid)
        slc_dict[uid]['sid'] = sid
    
    uids = list(slc_dict.keys())
    
    for uid in uids:
        t = find_in_annotation(annotation_dict=slc_dict[uid]['annotation'],
                               pattern='.//geolocationGridPoint/azimuthTime')
        swaths = t.keys()
        y = find_in_annotation(annotation_dict=slc_dict[uid]['annotation'], pattern='.//geolocationGridPoint/latitude')
        x = find_in_annotation(annotation_dict=slc_dict[uid]['annotation'], pattern='.//geolocationGridPoint/longitude')
        t_flat = np.asarray([datetime.fromisoformat(item).timestamp() for sublist in [t[swath] for swath in swaths]
                             for item in sublist])
        y_flat = np.asarray([float(item) for sublist in [y[swath] for swath in swaths] for item in sublist])
        x_flat = np.asarray([float(item) for sublist in [x[swath] for swath in swaths] for item in sublist])
        g = np.asarray([(x, y) for x, y in zip(x_flat, y_flat)])
        slc_dict[uid]['az_time'] = t_flat
        slc_dict[uid]['gridpts'] = g
    
    if len(uids) == 2:
        starts = [datetime.strptime(slc_dict[key]['sid'].start, '%Y%m%dT%H%M%S') for key in slc_dict.keys()]
        if starts[0] > starts[1]:
            az_time = np.concatenate([slc_dict[uids[1]]['az_time'], slc_dict[uids[0]]['az_time']])
            gridpts = np.concatenate([slc_dict[uids[1]]['gridpts'], slc_dict[uids[0]]['gridpts']])
        else:
            az_time = np.concatenate([slc_dict[key]['az_time'] for key in slc_dict.keys()])
            gridpts = np.concatenate([slc_dict[key]['gridpts'] for key in slc_dict.keys()])
    else:
        az_time = slc_dict[uids[0]]['az_time']
        gridpts = slc_dict[uids[0]]['gridpts']
    
    if slc_dict[uids[0]]['sid'].orbit == 'A':
        coord1 = lr
        coord2 = ul
    else:
        coord1 = ur
        coord2 = ll
    
    method = 'linear'
    res = [griddata(gridpts, az_time, coord1, method=method),
           griddata(gridpts, az_time, coord2, method=method)]
    
    min_start = min([datetime.strptime(slc_dict[uid]['sid'].start, '%Y%m%dT%H%M%S') for uid in uids])
    max_stop = max([datetime.strptime(slc_dict[uid]['sid'].stop, '%Y%m%dT%H%M%S') for uid in uids])
    res_t = []
    for i, r in enumerate(res):
        if np.isnan(r):
            if i == 0:
                res_t.append(min_start)
            else:
                res_t.append(max_stop)
        else:
            res_t.append(datetime.fromtimestamp(float(r)))
    
    start = datetime.strftime(res_t[0], '%Y%m%dT%H%M%S')
    stop = datetime.strftime(res_t[1], '%Y%m%dT%H%M%S')
    
    return start, stop


def create_data_mask(outname, datasets, extent, epsg, driver, creation_opt,
                     overviews, overview_resampling, dst_nodata, product_type, wbm=None):
    """
    Creation of the Data Mask image.
    
    Parameters
    ----------
    outname: str
        Full path to the output data mask file.
    datasets: list[dict]
        List of processed output files that match the source scenes and overlap with the current MGRS tile.
        An error will be thrown if not all datasets contain a key `datamask`.
        The function will return without an error if not all datasets contain a key `dm`.
    extent: dict
        Spatial extent of the MGRS tile, derived from a :class:`~spatialist.vector.Vector` object.
    epsg: int
        The coordinate reference system as an EPSG code.
    driver: str
        GDAL driver to use for raster file creation.
    creation_opt: list[str]
        GDAL creation options to use for raster file creation. Should match specified GDAL driver.
    overviews: list[int]
        Internal overview levels to be created for each raster file.
    overview_resampling: str
        Resampling method for overview levels.
    dst_nodata: int or str
        Nodata value to write to the output raster.
    product_type: str
        The type of ARD product that is being created. Either 'NRB' or 'ORB'.
    wbm: str or None
        Path to a water body mask file with the dimensions of an MGRS tile. Optional if `product_type='NRB', mandatory
        if `product_type='ORB'`.
    """
    measurement_keys = [x for x in datasets[0].keys() if re.search('[gs]-lin', x)]
    measurement = [scene[measurement_keys[0]] for scene in datasets]
    datamask = [scene['datamask'] for scene in datasets]
    ls = []
    for scene in datasets:
        if 'dm' in scene:
            ls.append(scene['dm'])
        else:
            return  # do not create a data mask if not all scenes have a layover-shadow mask
    
    dm_bands = ['not layover, nor shadow',
                'layover',
                'shadow',
                'ocean',
                'lakes',
                'rivers']
    
    if product_type == 'ORB':
        if wbm is None:
            raise RuntimeError('Water body mask is required for ORB products')
    
    tile_bounds = [extent['xmin'], extent['ymin'], extent['xmax'], extent['ymax']]
    
    vrt_ls = '/vsimem/' + os.path.dirname(outname) + 'ls.vrt'
    vrt_valid = '/vsimem/' + os.path.dirname(outname) + 'valid.vrt'
    vrt_measurement = '/vsimem/' + os.path.dirname(outname) + 'measurement.vrt'
    gdalbuildvrt(src=ls, dst=vrt_ls, outputBounds=tile_bounds, void=False)
    gdalbuildvrt(src=datamask, dst=vrt_valid, outputBounds=tile_bounds, void=False)
    gdalbuildvrt(src=measurement, dst=vrt_measurement, outputBounds=tile_bounds, void=False)
    
    with Raster(vrt_ls) as ras_ls:
        with bbox(extent, crs=epsg) as tile_vec:
            ras_ls_res = ras_ls.res
            rows = ras_ls.rows
            cols = ras_ls.cols
            geotrans = ras_ls.raster.GetGeoTransform()
            proj = ras_ls.raster.GetProjection()
            arr_dm = ras_ls.array()
            
            # Get Water Body Mask
            if wbm is not None:
                with Raster(wbm) as ras_wbm:
                    ras_wbm_cols = ras_wbm.cols
                    cols_ratio = ras_wbm_cols / cols
                    if cols_ratio > 1:
                        # create low resolution VRT
                        res = int(ras_ls_res[0])
                        wbm_lowres = wbm.replace('.tif', f'_{res}m.vrt')
                        options = {'xRes': res, 'yRes': res,
                                   'resampleAlg': 'mode'}
                        gdalbuildvrt(src=wbm, dst=wbm_lowres, **options)
                        with Raster(wbm_lowres) as ras_wbm_lowres:
                            arr_wbm = ras_wbm_lowres.array()
                    else:
                        arr_wbm = ras_wbm.array()
            else:
                del dm_bands[3:]
            
            # Extend the shadow class of the data mask with nodata values
            # from backscatter data and create final array
            with Raster(vrt_valid)[tile_vec] as ras_valid:
                with Raster(vrt_measurement)[tile_vec] as ras_measurement:
                    arr_valid = ras_valid.array()
                    arr_measurement = ras_measurement.array()
                    
                    arr_dm = np.nan_to_num(arr_dm)
                    arr_dm = np.where(((arr_valid == 1) & (np.isnan(arr_measurement))),
                                      2, arr_dm)
                    arr_dm[np.isnan(arr_valid)] = dst_nodata
                    del arr_measurement
                    del arr_valid
        
        outname_tmp = '/vsimem/' + os.path.basename(outname) + '.vrt'
        gdriver = gdal.GetDriverByName('GTiff')
        ds_tmp = gdriver.Create(outname_tmp, rows, cols, len(dm_bands), gdal.GDT_Byte,
                                options=['ALPHA=UNSPECIFIED', 'PHOTOMETRIC=MINISWHITE'])
        gdriver = None
        ds_tmp.SetGeoTransform(geotrans)
        ds_tmp.SetProjection(proj)
        
        for i, name in enumerate(dm_bands):
            band = ds_tmp.GetRasterBand(i + 1)
            
            # not layover, nor shadow
            if i == 0:
                arr = arr_dm == i
            # layover | shadow
            # source value 3: layover and shadow
            elif i in [1, 2]:
                arr = (arr_dm == i) | (arr_dm == 3)
            # ocean
            elif i == 3:
                arr = arr_wbm == 1
            # lakes
            elif i == 4:
                arr = arr_wbm == 2
            # rivers
            elif i == 5:
                arr = arr_wbm == 3
            else:
                raise ValueError(f'unknown array value: {i}')
            
            arr = arr.astype('uint8')
            band.WriteArray(arr)
            band.SetNoDataValue(dst_nodata)
            band.SetDescription(name)
            band.FlushCache()
            band = None
            del arr
        
        ds_tmp.SetMetadataItem('TIFFTAG_DATETIME', strftime('%Y:%m:%d %H:%M:%S', gmtime()))
        ds_tmp.BuildOverviews(overview_resampling, overviews)
        outDataset_cog = gdal.GetDriverByName(driver).CreateCopy(outname, ds_tmp, strict=1, options=creation_opt)
        outDataset_cog = None
        ds_tmp = None
        tile_vec = None


def create_acq_id_image(outname, ref_tif, datasets, src_ids, extent,
                        epsg, driver, creation_opt, overviews, dst_nodata):
    """
    Creation of the Acquisition ID image.

    Parameters
    ----------
    outname: str
        Full path to the output data mask file.
    ref_tif: str
        Full path to any GeoTIFF file of the ARD product.
    datasets: list[dict]
        List of processed output files that match the source SLC scenes and overlap with the current MGRS tile.
    src_ids: list[pyroSAR.drivers.ID]
        List of :class:`~pyroSAR.drivers.ID` objects of all source SLC scenes that overlap with the current MGRS tile.
    extent: dict
        Spatial extent of the MGRS tile, derived from a :class:`~spatialist.vector.Vector` object.
    epsg: int
        The CRS used for the ARD product; provided as an EPSG code.
    driver: str
        GDAL driver to use for raster file creation.
    creation_opt: list[str]
        GDAL creation options to use for raster file creation. Should match specified GDAL driver.
    overviews: list[int]
        Internal overview levels to be created for each raster file.
    dst_nodata: int or str
        Nodata value to write to the output raster.
    """
    src_scenes = [sid.scene for sid in src_ids]
    # If there are two source scenes, make sure that the order of acquisitions in all lists is correct!
    if len(src_scenes) > 1:
        if not len(src_scenes) == 2 and len(datasets) == 2:
            raise RuntimeError('expected lists `src_scenes` and `valid_mask_list` to be of length 2; length is '
                               '{} and {} respectively'.format(len(src_scenes), len(datasets)))
        starts_src = [datetime.strptime(identify(f).start, '%Y%m%dT%H%M%S') for f in src_scenes]
        start_valid = [re.search('[0-9]{8}T[0-9]{6}', os.path.basename(x)).group() for x in src_scenes]
        start_valid = [datetime.strptime(x, '%Y%m%dT%H%M%S') for x in start_valid]
        if starts_src[0] > starts_src[1]:
            src_scenes.reverse()
            starts_src.reverse()
        if start_valid[0] != starts_src[0]:
            datasets.reverse()
        if start_valid[0] != starts_src[0]:
            raise RuntimeError('failed to match order of lists `src_scenes` and `valid_mask_list`')
    
    tile_bounds = [extent['xmin'], extent['ymin'], extent['xmax'], extent['ymax']]
    
    arr_list = []
    for dataset in datasets:
        vrt_valid = '/vsimem/' + os.path.dirname(outname) + 'mosaic.vrt'
        gdalbuildvrt(src=dataset['datamask'], dst=vrt_valid, outputBounds=tile_bounds, void=False)
        with bbox(extent, crs=epsg) as tile_vec:
            with Raster(vrt_valid)[tile_vec] as vrt_ras:
                vrt_arr = vrt_ras.array()
                arr_list.append(vrt_arr)
                del vrt_arr
            tile_vec = None
    
    src_scenes_clean = [os.path.basename(src).replace('.zip', '').replace('.SAFE', '') for src in src_scenes]
    tag = '{{"{src1}": 1}}'.format(src1=src_scenes_clean[0])
    out_arr = np.full(arr_list[0].shape, dst_nodata)
    out_arr[arr_list[0] == 1] = 1
    if len(arr_list) == 2:
        out_arr[arr_list[1] == 1] = 2
        tag = '{{"{src1}": 1, "{src2}": 2}}'.format(src1=src_scenes_clean[0], src2=src_scenes_clean[1])
    
    creation_opt.append('TIFFTAG_IMAGEDESCRIPTION={}'.format(tag))
    with Raster(ref_tif) as ref_ras:
        ref_ras.write(outname, format=driver, array=out_arr.astype('uint8'), nodata=dst_nodata, overwrite=True,
                      overviews=overviews, options=creation_opt)


def wind_normalization(src, dst_wm, dst_wn, measurement, gapfill, bounds, epsg, driver, creation_opt, dst_nodata,
                       multithread, resolution=915):
    """
    Create wind normalization layers. A wind model annotation layer is created and optionally
    a wind normalization VRT.
    
    Parameters
    ----------
    src: list[str]
        A list of OCN products as prepared by :func:`s1ard.ocn.extract`
    dst_wm: str
        The name of the wind model layer in the ARD product
    dst_wn: str or None
        The name of the wind normalization VRT. If None, no VRT will be created.
        Requires `measurement` to point to a file.
    measurement: str or None
        The name of the measurement file used for wind normalization in `dst_wn`.
        If None, no wind normalization VRT will be created.
    gapfill: bool
        Perform additional gap filling (:func:`s1ard.ocn.gapfill`)?
        This is recommended if the Level-1 source product of `measurement` is GRD
        in which case gaps are introduced between subsequently acquired scenes.
    bounds: list[float]
        the bounds of the MGRS tile
    epsg: int
        The EPSG code of the MGRS tile
    driver: str
        GDAL driver to use for raster file creation.
    creation_opt: list[str]
        GDAL creation options to use for raster file creation. Should match specified GDAL driver.
    dst_nodata: float
        Nodata value to write to the output raster.
    multithread: bool
        Should `gdalwarp` use multithreading?
    resolution: int, optional
        The target pixel resolution in meters. 915 is chosen as default because it is closest
        to the OCN product resolution (1000) and still fits into the MGRS bounds
        (``109800 % 915 == 0``).
    
    Returns
    -------
    
    """
    if len(src) > 1:
        cmod_mosaic = tempfile.NamedTemporaryFile(suffix='.tif').name
        gdalwarp(src=src, dst=cmod_mosaic)
        if gapfill:
            cmod_geo = tempfile.NamedTemporaryFile(suffix='.tif').name
            ocn.gapfill(src=cmod_mosaic, dst=cmod_geo, md=2, si=1)
        else:
            cmod_geo = cmod_mosaic
    else:
        cmod_geo = src[0]
    
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
    
    if dst_wn is not None and measurement is not None:
        if not os.path.isfile(dst_wn):
            log.info(f'creating {dst_wn}')
            with Raster(measurement) as ras:
                xres, yres = ras.res
            create_vrt(src=[measurement, dst_wm], dst=dst_wn, fun='div',
                       options={'xRes': xres, 'yRes': yres}, relpaths=True)
