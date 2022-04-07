import os
import re
import shutil
import tempfile
from datetime import datetime, timezone
from lxml import etree
import numpy as np
from spatialist import Raster, Vector, vectorize, boundary, bbox, intersect, rasterize
from spatialist.auxil import gdalwarp, gdalbuildvrt
from pyroSAR import identify_many
from pyroSAR.ancillary import find_datasets
import S1_NRB
from S1_NRB.metadata import extract, xmlparser, stacparser
from S1_NRB.metadata.mapping import ITEM_MAP
import S1_NRB.ancillary as ancil


def format(config, scenes, datadir, outdir, tile, extent, epsg, wbm=None,
           multithread=True, compress=None, overviews=None):
    """
    Finalizes the generation of Sentinel-1 NRB products after processing steps via `pyroSAR.snap.util.geocode` and
    `pyroSAR.snap.util.noise_power` have been finished. This includes the following:
    - Creating all measurement and annotation datasets in Cloud Optimized GeoTIFF (COG) format
    - Creating all annotation datasets in Virtual Raster Tile (VRT) format
    - Applying the NRB product directory structure & naming convention
    - Generating metadata in XML and JSON formats for the NRB product as well as source SLC datasets

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

    Returns
    -------
    None
    """
    if compress is None:
        compress = 'LERC_ZSTD'
    if overviews is None:
        overviews = [2, 4, 9, 18, 36]
    ovr_resampling = 'AVERAGE'
    driver = 'COG'
    blocksize = 512
    src_nodata_snap = 0
    dst_nodata = 'nan'
    dst_nodata_mask = 255
    
    src_ids, snap_datasets, snap_datamasks = get_datasets(scenes=scenes, datadir=datadir,
                                                          tile=tile, extent=extent, epsg=epsg)
    nrb_start, nrb_stop = ancil.calc_product_start_stop(src_ids=src_ids, extent=extent, epsg=epsg)
    meta = {'mission': src_ids[0].sensor,
            'mode': src_ids[0].meta['acquisition_mode'],
            'polarization': {"['HH']": 'SH',
                             "['VV']": 'SV',
                             "['HH', 'HV']": 'DH',
                             "['VV', 'VH']": 'DV'}[str(src_ids[0].polarizations)],
            'start': nrb_start,
            'orbitnumber': src_ids[0].meta['orbitNumbers_abs']['start'],
            'datatake': hex(src_ids[0].meta['frameNumber']).replace('x', '').upper(),
            'tile': tile,
            'id': 'ABCD'}
    meta_lower = dict((k, v.lower() if not isinstance(v, int) else v) for k, v in meta.items())
    skeleton_dir = '{mission}_{mode}_NRB__1S{polarization}_{start}_{orbitnumber:06}_{datatake}_{tile}_{id}'
    skeleton_files = '{mission}-{mode}-nrb-{start}-{orbitnumber:06}-{datatake}-{tile}-{suffix}.tif'
    
    nrb_dir = os.path.join(outdir, skeleton_dir.format(**meta))
    os.makedirs(nrb_dir, exist_ok=True)
    
    # prepare raster write options; https://gdal.org/drivers/raster/cog.html
    write_options_base = ['BLOCKSIZE={}'.format(blocksize), 'OVERVIEW_RESAMPLING={}'.format(ovr_resampling)]
    write_options = dict()
    for key in ITEM_MAP:
        write_options[key] = write_options_base.copy()
        if compress is not None:
            entry = 'COMPRESS={}'.format(compress)
            write_options[key].append(entry)
            if compress.startswith('LERC'):
                entry = 'MAX_Z_ERROR={:f}'.format(ITEM_MAP[key]['z_error'])
                write_options[key].append(entry)
    
    # create raster files: linear gamma0 backscatter (-[vh|vv|hh|hv]-g-lin.tif), ellipsoidal incident angle (-ei.tif),
    # gamma-to-sigma ratio (-gs.tif), local contributing area (-lc.tif), local incident angle (-li.tif),
    # noise power images (-np-[vh|vv|hh|hv].tif)
    nrb_tifs = []
    pattern = '|'.join(ITEM_MAP.keys())
    for i, ds in enumerate(snap_datasets):
        if isinstance(ds, str):
            match = re.search(pattern, ds)
            if match is not None:
                key = match.group()
            else:
                continue
        else:
            match = [re.search(pattern, x) for x in ds]
            keys = [x if x is None else x.group() for x in match]
            if len(list(set(keys))) != 1:
                raise RuntimeError('file mismatch:\n{}'.format('\n'.join(ds)))
            if None in keys:
                continue
            key = keys[0]
        
        if key == 'layoverShadowMask':
            # the data mask raster (-dm.tif) will be created later on in the processing workflow
            continue
        
        meta_lower['suffix'] = ITEM_MAP[key]['suffix']
        outname_base = skeleton_files.format(**meta_lower)
        if re.search('_gamma0', key):
            subdir = 'measurement'
        else:
            subdir = 'annotation'
        outname = os.path.join(nrb_dir, subdir, outname_base)
        
        if not os.path.isfile(outname):
            os.makedirs(os.path.dirname(outname), exist_ok=True)
            print(outname)
            bounds = [extent['xmin'], extent['ymin'], extent['xmax'], extent['ymax']]
            
            if isinstance(ds, tuple):
                with Raster(list(ds), list_separate=False) as ras:
                    source = ras.filename
            elif isinstance(ds, str):
                source = tempfile.NamedTemporaryFile(suffix='.vrt').name
                gdalbuildvrt(ds, source)
            else:
                raise TypeError('type {} is not supported: {}'.format(type(ds), ds))
            
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
            
            gdalwarp(source, outname,
                     options={'format': driver, 'outputBounds': bounds, 'srcNodata': src_nodata_snap,
                              'dstNodata': dst_nodata, 'multithread': multithread,
                              'creationOptions': write_options[key]})
            nrb_tifs.append(outname)
    
    # determine processing timestamp, generate unique ID and rename NRB directory accordingly
    proc_time = datetime.now(timezone.utc)
    t = proc_time.isoformat().encode()
    product_id = ancil.generate_unique_id(encoded_str=t)
    nrb_dir_new = nrb_dir.replace('ABCD', product_id)
    nrb_tifs = [tif.replace(nrb_dir, nrb_dir_new) for tif in nrb_tifs]
    os.rename(nrb_dir, nrb_dir_new)
    nrb_dir = nrb_dir_new
    
    # reformat `snap_datasets` to a flattened list if necessary
    if type(snap_datasets[0]) == tuple:
        snap_datasets = [item for tup in snap_datasets for item in tup]
    
    # define a reference raster from the annotation datasets and list all gamma0 backscatter measurement rasters
    ref_tif_suffix = '-lc.tif'
    ref_tif = [tif for tif in nrb_tifs if re.search('{}$'.format(ref_tif_suffix), tif)][0]
    measure_tifs = [tif for tif in nrb_tifs if re.search('[hv]{2}-g-lin.tif$', tif)]
    
    # create data mask raster (-dm.tif)
    if wbm is not None:
        if not config['dem_type'] == 'GETASSE30' and not os.path.isfile(wbm):
            raise FileNotFoundError('External water body mask could not be found: {}'.format(wbm))
    
    dm_path = ref_tif.replace(ref_tif_suffix, '-dm.tif')
    ancil.create_data_mask(outname=dm_path, snap_datamasks=snap_datamasks, snap_datasets=snap_datasets,
                           extent=extent, epsg=epsg, driver=driver, creation_opt=write_options['layoverShadowMask'],
                           overviews=overviews, overview_resampling=ovr_resampling, wbm=wbm, dst_nodata=dst_nodata_mask)
    nrb_tifs.append(dm_path)
    
    # create acquisition ID image raster (-id.tif)
    id_path = ref_tif.replace(ref_tif_suffix, '-id.tif')
    ancil.create_acq_id_image(outname=id_path, ref_tif=ref_tif, snap_datamasks=snap_datamasks, src_ids=src_ids,
                              extent=extent, epsg=epsg, driver=driver, creation_opt=write_options['acquisitionImage'],
                              overviews=overviews, dst_nodata=dst_nodata_mask)
    nrb_tifs.append(id_path)
    
    # create color composite VRT (-cc-g-lin.vrt)
    if meta['polarization'] in ['DH', 'DV'] and len(measure_tifs) == 2:
        cc_path = re.sub('[hv]{2}', 'cc', measure_tifs[0]).replace('.tif', '.vrt')
        ancil.create_rgb_vrt(outname=cc_path, infiles=measure_tifs, overviews=overviews,
                             overview_resampling=ovr_resampling)
    
    # create log-scaled gamma nought VRTs (-[vh|vv|hh|hv]-g-log.vrt)
    for item in measure_tifs:
        gamma0_rtc_log = item.replace('lin.tif', 'log.vrt')
        if not os.path.isfile(gamma0_rtc_log):
            print(gamma0_rtc_log)
            ancil.create_vrt(src=item, dst=gamma0_rtc_log, fun='log10', scale=10,
                             options={'VRTNodata': 'NaN'}, overviews=overviews, overview_resampling=ovr_resampling)
    
    # create sigma nought RTC VRTs (-[vh|vv|hh|hv]-s-[lin|log].vrt)
    for item in measure_tifs:
        sigma0_rtc_lin = item.replace('g-lin.tif', 's-lin.vrt')
        sigma0_rtc_log = item.replace('g-lin.tif', 's-log.vrt')
        
        if not os.path.isfile(sigma0_rtc_lin):
            print(sigma0_rtc_lin)
            ancil.create_vrt(src=[item, ref_tif], dst=sigma0_rtc_lin, fun='mul', relpaths=True,
                             options={'VRTNodata': 'NaN'}, overviews=overviews, overview_resampling=ovr_resampling)
        
        if not os.path.isfile(sigma0_rtc_log):
            print(sigma0_rtc_log)
            ancil.create_vrt(src=sigma0_rtc_lin, dst=sigma0_rtc_log, fun='log10', scale=10,
                             options={'VRTNodata': 'NaN'}, overviews=overviews, overview_resampling=ovr_resampling)
    
    # copy support files
    schema_dir = os.path.join(S1_NRB.__path__[0], 'validation', 'schemas')
    schemas = [os.path.join(schema_dir, schema) for schema in os.listdir(schema_dir)]
    support_dir = os.path.join(nrb_dir, 'support')
    os.makedirs(support_dir, exist_ok=True)
    for schema in schemas:
        shutil.copy(schema, support_dir)
    
    # create metadata files in XML and (STAC) JSON formats
    meta = extract.meta_dict(config=config, target=nrb_dir, src_ids=src_ids, snap_datasets=snap_datasets,
                             proc_time=proc_time, start=nrb_start, stop=nrb_stop, compression=compress)
    xmlparser.main(meta=meta, target=nrb_dir, tifs=nrb_tifs)
    stacparser.main(meta=meta, target=nrb_dir, tifs=nrb_tifs)


def get_datasets(scenes, datadir, tile, extent, epsg):
    """
    Identifies all source SLC scenes, finds matching output files processed with `pyroSAR.snap.util.geocode` in
    `datadir` and filters both lists depending on the actual overlap of each SLC footprint with the current MGRS tile
    geometry.

    Parameters
    ----------
    scenes: list[str]
        List of scenes to process. Either an individual scene or multiple, matching scenes (consecutive acquisitions).
    datadir: str
        The directory containing the datasets processed from the source scenes using pyroSAR.
    tile: str
        ID of an MGRS tile.
    extent: dict
        Spatial extent of the MGRS tile, derived from a `spatialist.vector.Vector` object.
    epsg: int
        The CRS used for the NRB product; provided as an EPSG code.

    Returns
    -------
    ids: list[ID]
        List of `pyroSAR.driver.ID` objects of all source SLC scenes that overlap with the current MGRS tile.
    datasets: list[str] or list[list[str]]
        List of output files processed with `pyroSAR.snap.util.geocode` that match each ID object of `ids`.
        The format of `snap_datasets` is a list of strings if only a single ID object is stored in `ids`, else it is
        a list of lists.
    datamasks: list[str]
        List of raster datamask files covering the footprint of each source SLC scene that overlaps with the current
        MGRS tile.
    """
    ids = identify_many(scenes)
    datasets = []
    for _id in ids:
        scene_base = os.path.splitext(os.path.basename(_id.scene))[0]
        scene_dir = os.path.join(datadir, scene_base, str(epsg))
        datasets.append(find_datasets(directory=scene_dir))
    if len(datasets) == 0:
        raise RuntimeError("No pyroSAR datasets were found in the directory '{}'".format(datadir))
    
    pattern = '[VH]{2}_gamma0-rtc'
    i = 0
    datamasks = []
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
                    datamasks.append(snap_dm_ras)
                    i += 1
                    inter.close()
    if len(ids) == 0:
        raise RuntimeError('None of the scenes overlap with the current tile {tile_id}: '
                           '\n{scenes}'.format(tile_id=tile, scenes=scenes))
    
    if len(datasets) > 1:
        datasets = list(zip(*datasets))
    else:
        datasets = datasets[0]
    
    return ids, datasets, datamasks
