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
from spatialist import Raster, Vector, vectorize, boundary, bbox, intersect, rasterize
from spatialist.auxil import gdalwarp, gdalbuildvrt
from spatialist.ancillary import finder
from pyroSAR import identify, identify_many
from pyroSAR.ancillary import find_datasets
import S1_NRB
from S1_NRB.metadata import extract, xml, stac
from S1_NRB.metadata.mapping import ITEM_MAP
from S1_NRB.ancillary import generate_unique_id
from S1_NRB.metadata.extract import etree_from_sid, find_in_annotation


def format(config, scenes, datadir, outdir, tile, extent, epsg, wbm=None,
           multithread=True, compress=None, overviews=None):
    """
    Finalizes the generation of Sentinel-1 NRB products after processing steps via :func:`pyroSAR.snap.util.geocode` and
    :func:`pyroSAR.snap.util.noise_power` have finished. This includes the following:
    - Creating all measurement and annotation datasets in Cloud Optimized GeoTIFF (COG) format
    - Creating additional annotation datasets in Virtual Raster Tile (VRT) format
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
        Spatial extent of the MGRS tile, derived from a :class:`~spatialist.vector.Vector` object.
    epsg: int
        The CRS used for the NRB product; provided as an EPSG code.
    wbm: str, optional
        Path to a water body mask file with the dimensions of an MGRS tile.
    multithread: bool, optional
        Should `gdalwarp` use multithreading? Default is True. The number of threads used, can be adjusted in the
        `config.ini` file with the parameter `gdal_threads`.
    compress: str, optional
        Compression algorithm to use. See https://gdal.org/drivers/raster/gtiff.html#creation-options for options.
        Defaults to 'LERC_DEFLATE'.
    overviews: list[int], optional
        Internal overview levels to be created for each GeoTIFF file. Defaults to [2, 4, 9, 18, 36]

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
    src_nodata = 0
    dst_nodata_float = -9999.0
    dst_nodata_byte = 255
    vrt_nodata = 'nan'  # was found necessary for proper calculation of statistics in QGIS
    
    # determine processing timestamp and generate unique ID
    start_time = time.time()
    proc_time = datetime.now(timezone.utc)
    t = proc_time.isoformat().encode()
    product_id = generate_unique_id(encoded_str=t)
    
    src_ids, snap_datasets, snap_datamasks = get_datasets(scenes=scenes, datadir=datadir,
                                                          tile=tile, extent=extent, epsg=epsg)
    nrb_start, nrb_stop = calc_product_start_stop(src_ids=src_ids, extent=extent, epsg=epsg)
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
            'id': product_id}
    meta_lower = dict((k, v.lower() if not isinstance(v, int) else v) for k, v in meta.items())
    skeleton_dir = '{mission}_{mode}_NRB__1S{polarization}_{start}_{orbitnumber:06}_{datatake:0>6}_{tile}_{id}'
    skeleton_files = '{mission}-{mode}-nrb-{start}-{orbitnumber:06}-{datatake:0>6}-{tile}-{suffix}.tif'
    
    modify_existing = False  # modify existing products so that only missing files are re-created
    nrb_base = skeleton_dir.format(**meta)
    existing = finder(outdir, [nrb_base.replace(product_id, '*')], foldermode=2)
    if len(existing) > 0:
        if not modify_existing:
            return 'Already processed - Skip!'
        else:
            nrb_dir = existing[0]
    else:
        nrb_dir = os.path.join(outdir, nrb_base)
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
            
            ras = None
            if isinstance(ds, tuple):
                ras = Raster(list(ds), list_separate=False)
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
                     options={'format': driver, 'outputBounds': bounds, 'srcNodata': src_nodata,
                              'dstNodata': dst_nodata_float, 'multithread': multithread,
                              'creationOptions': write_options[key]})
            if ras is not None:
                ras.close()
        nrb_tifs.append(outname)
    
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
    if not os.path.isfile(dm_path):
        create_data_mask(outname=dm_path, snap_datamasks=snap_datamasks,
                         snap_datasets=snap_datasets, extent=extent, epsg=epsg,
                         driver=driver, creation_opt=write_options['layoverShadowMask'],
                         overviews=overviews, overview_resampling=ovr_resampling,
                         wbm=wbm, dst_nodata=dst_nodata_byte)
    nrb_tifs.append(dm_path)
    
    # create acquisition ID image raster (-id.tif)
    id_path = ref_tif.replace(ref_tif_suffix, '-id.tif')
    if not os.path.isfile(id_path):
        create_acq_id_image(outname=id_path, ref_tif=ref_tif,
                            snap_datamasks=snap_datamasks, src_ids=src_ids,
                            extent=extent, epsg=epsg, driver=driver,
                            creation_opt=write_options['acquisitionImage'],
                            overviews=overviews, dst_nodata=dst_nodata_byte)
    nrb_tifs.append(id_path)
    
    # create color composite VRT (-cc-g-lin.vrt)
    if meta['polarization'] in ['DH', 'DV'] and len(measure_tifs) == 2:
        cc_path = re.sub('[hv]{2}', 'cc', measure_tifs[0]).replace('.tif', '.vrt')
        if not os.path.isfile(cc_path):
            create_rgb_vrt(outname=cc_path, infiles=measure_tifs,
                           overviews=overviews, overview_resampling=ovr_resampling)
    
    vrt_options = {'VRTNodata': vrt_nodata}
    
    # create log-scaled gamma nought VRTs (-[vh|vv|hh|hv]-g-log.vrt)
    fun = 'dB'
    args = {'fact': 10}
    scale = None
    for item in measure_tifs:
        gamma0_rtc_log = item.replace('lin.tif', 'log.vrt')
        if not os.path.isfile(gamma0_rtc_log):
            print(gamma0_rtc_log)
            create_vrt(src=item, dst=gamma0_rtc_log, fun=fun, scale=scale,
                       args=args, options=vrt_options, overviews=overviews,
                       overview_resampling=ovr_resampling)
    
    # create sigma nought RTC VRTs (-[vh|vv|hh|hv]-s-[lin|log].vrt)
    gs_path = ref_tif.replace(ref_tif_suffix, '-gs.tif')
    for item in measure_tifs:
        sigma0_rtc_lin = item.replace('g-lin.tif', 's-lin.vrt')
        sigma0_rtc_log = item.replace('g-lin.tif', 's-log.vrt')
        
        if not os.path.isfile(sigma0_rtc_lin):
            print(sigma0_rtc_lin)
            create_vrt(src=[item, gs_path], dst=sigma0_rtc_lin, fun='mul',
                       relpaths=True, options=vrt_options, overviews=overviews,
                       overview_resampling=ovr_resampling)
        
        if not os.path.isfile(sigma0_rtc_log):
            print(sigma0_rtc_log)
            create_vrt(src=sigma0_rtc_lin, dst=sigma0_rtc_log, fun=fun,
                       scale=scale, options=vrt_options, overviews=overviews,
                       overview_resampling=ovr_resampling, args=args)
    
    # copy support files
    schema_dir = os.path.join(S1_NRB.__path__[0], 'validation', 'schemas')
    schemas = [os.path.join(schema_dir, schema) for schema in os.listdir(schema_dir)]
    support_dir = os.path.join(nrb_dir, 'support')
    os.makedirs(support_dir, exist_ok=True)
    for schema in schemas:
        shutil.copy(schema, support_dir)
    
    # create metadata files in XML and (STAC) JSON formats
    start = datetime.strptime(nrb_start, '%Y%m%dT%H%M%S')
    stop = datetime.strptime(nrb_stop, '%Y%m%dT%H%M%S')
    meta = extract.meta_dict(config=config, target=nrb_dir, src_ids=src_ids, snap_datasets=snap_datasets,
                             proc_time=proc_time, start=start, stop=stop, compression=compress)
    xml.parse(meta=meta, target=nrb_dir, tifs=nrb_tifs, exist_ok=True)
    stac.parse(meta=meta, target=nrb_dir, tifs=nrb_tifs, exist_ok=True)
    return str(round((time.time() - start_time), 2))


def get_datasets(scenes, datadir, tile, extent, epsg):
    """
    Identifies all source SLC scenes, finds matching output files processed with :func:`pyroSAR.snap.util.geocode` in
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
        Spatial extent of the MGRS tile, derived from a :class:`~spatialist.vector.Vector` object.
    epsg: int
        The coordinate reference system as an EPSG code.

    Returns
    -------
    ids: list[:class:`pyroSAR.drivers.ID`]
        List of :class:`~pyroSAR.drivers.ID` objects of all source SLC scenes that overlap with the current MGRS tile.
    datasets: list[str] or list[list[str]]
        List of output files processed with :func:`pyroSAR.snap.util.geocode` that match each
        :class:`~pyroSAR.drivers.ID` object of `ids`. The format is a list of strings if only a single object is stored
        in `ids`, else it is a list of lists.
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


def create_vrt(src, dst, fun, relpaths=False, scale=None, offset=None, args=None,
               options=None, overviews=None, overview_resampling=None):
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
    relpaths: bool, optional
        Should all `SourceFilename` XML elements with attribute `@relativeToVRT="0"` be updated to be paths relative to
        the output VRT file? Default is False.
    scale: int, optional
         The scale that should be applied when computing “real” pixel values from scaled pixel values on a raster band.
         Will be ignored if `fun='decibel'`.
    offset: float, optional
        The offset that should be applied when computing “real” pixel values from scaled pixel values on a raster band.
        Will be ignored if `fun='decibel'`.
    args: dict, optional
        arguments for `fun` passed as `PixelFunctionArguments`. Requires GDAL>=3.5 to be read.
    options: dict, optional
        Additional parameters passed to `gdal.BuildVRT`.
    overviews: list[int], optional
        Internal overview levels to be created for each raster file.
    overview_resampling: str, optional
        Resampling method for overview levels.

    Examples
    --------
    linear backscatter as input:

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
    gdalbuildvrt(src=src, dst=dst, options=options)
    tree = etree.parse(dst)
    root = tree.getroot()
    band = tree.find('VRTRasterBand')
    band.attrib['subClass'] = 'VRTDerivedRasterBand'
    
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
    print(outname)
    
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
    gdalbuildvrt(src=infiles, dst=outname,
                 options={'separate': True})
    
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
    Calculates the start and stop times of the NRB product.
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
    start: str
        Start time of the NRB product formatted as `%Y%m%dT%H%M%S` in UTC.
    stop: str
        Stop time of the NRB product formatted as `%Y%m%dT%H%M%S` in UTC.
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
        slc_dict[uid] = etree_from_sid(sid)
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


def create_data_mask(outname, snap_datamasks, snap_datasets, extent, epsg, driver,
                     creation_opt, overviews, overview_resampling, dst_nodata, wbm=None):
    """
    Creation of the Data Mask image.

    Parameters
    ----------
    outname: str
        Full path to the output data mask file.
    snap_datamasks: list[str]
        List of raster datamask files covering the footprint of each source SLC scene that overlaps with the current
        MGRS tile.
    snap_datasets: list[str]
        List of output files processed with :func:`pyroSAR.snap.util.geocode` that match the source SLC scenes and
        overlap with the current MGRS tile.
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
    wbm: str, optional
        Path to a water body mask file with the dimensions of an MGRS tile.
    """
    print(outname)
    pols = [pol for pol in set([re.search('[VH]{2}', os.path.basename(x)).group() for x in snap_datasets if
                                re.search('[VH]{2}', os.path.basename(x)) is not None])]
    pattern = pols[0] + '_gamma0-rtc'
    snap_gamma0 = [x for x in snap_datasets if re.search(pattern, os.path.basename(x))]
    snap_ls_mask = [x for x in snap_datasets if re.search('layoverShadowMask', os.path.basename(x))]
    
    dm_bands = {1: {'arr_val': 0,
                    'name': 'not layover, nor shadow'},
                2: {'arr_val': 1,
                    'name': 'layover'},
                3: {'arr_val': 2,
                    'name': 'shadow'},
                4: {'arr_val': 4,
                    'name': 'ocean water'}}
    
    tile_bounds = [extent['xmin'], extent['ymin'], extent['xmax'], extent['ymax']]
    
    vrt_snap_ls = '/vsimem/' + os.path.dirname(outname) + 'snap_ls.vrt'
    vrt_snap_valid = '/vsimem/' + os.path.dirname(outname) + 'snap_valid.vrt'
    vrt_snap_gamma0 = '/vsimem/' + os.path.dirname(outname) + 'snap_gamma0.vrt'
    gdalbuildvrt(snap_ls_mask, vrt_snap_ls, options={'outputBounds': tile_bounds}, void=False)
    gdalbuildvrt(snap_datamasks, vrt_snap_valid, options={'outputBounds': tile_bounds}, void=False)
    gdalbuildvrt(snap_gamma0, vrt_snap_gamma0, options={'outputBounds': tile_bounds}, void=False)
    
    with Raster(vrt_snap_ls) as ras_snap_ls:
        with bbox(extent, crs=epsg) as tile_vec:
            rows = ras_snap_ls.rows
            cols = ras_snap_ls.cols
            geotrans = ras_snap_ls.raster.GetGeoTransform()
            proj = ras_snap_ls.raster.GetProjection()
            arr_snap_dm = ras_snap_ls.array()
            
            # Add Water Body Mask
            if wbm is not None:
                with Raster(wbm) as ras_wbm:
                    arr_wbm = ras_wbm.array()
                    out_arr = np.where((arr_wbm == 1), 4, arr_snap_dm)
                    del arr_wbm
            else:
                out_arr = arr_snap_dm
                dm_bands.pop(4)
            del arr_snap_dm
            
            # Extend the shadow class of the data mask with nodata values from backscatter data and create final array
            with Raster(vrt_snap_valid)[tile_vec] as ras_snap_valid:
                with Raster(vrt_snap_gamma0)[tile_vec] as ras_snap_gamma0:
                    arr_snap_valid = ras_snap_valid.array()
                    arr_snap_gamma0 = ras_snap_gamma0.array()
                    
                    out_arr = np.nan_to_num(out_arr)
                    out_arr = np.where(((arr_snap_valid == 1) & (np.isnan(arr_snap_gamma0)) & (out_arr != 4)), 2,
                                       out_arr)
                    out_arr[np.isnan(arr_snap_valid)] = dst_nodata
                    del arr_snap_gamma0
                    del arr_snap_valid
        
        outname_tmp = '/vsimem/' + os.path.basename(outname) + '.vrt'
        gdriver = gdal.GetDriverByName('GTiff')
        ds_tmp = gdriver.Create(outname_tmp, rows, cols, len(dm_bands.keys()), gdal.GDT_Byte,
                                options=['ALPHA=UNSPECIFIED', 'PHOTOMETRIC=MINISWHITE'])
        gdriver = None
        ds_tmp.SetGeoTransform(geotrans)
        ds_tmp.SetProjection(proj)
        
        for k, v in dm_bands.items():
            band = ds_tmp.GetRasterBand(k)
            arr_val = v['arr_val']
            b_name = v['name']
            
            arr = np.full((rows, cols), 0)
            arr[out_arr == dst_nodata] = dst_nodata
            if arr_val == 0:
                arr[out_arr == 0] = 1
            elif arr_val in [1, 2]:
                arr[(out_arr == arr_val) | (out_arr == 3)] = 1
            elif arr_val == 4:
                arr[out_arr == 4] = 1
            
            arr = arr.astype('uint8')
            band.WriteArray(arr)
            band.SetNoDataValue(dst_nodata)
            band.SetDescription(b_name)
            band.FlushCache()
            band = None
            del arr
        
        ds_tmp.SetMetadataItem('TIFFTAG_DATETIME', strftime('%Y:%m:%d %H:%M:%S', gmtime()))
        ds_tmp.BuildOverviews(overview_resampling, overviews)
        outDataset_cog = gdal.GetDriverByName(driver).CreateCopy(outname, ds_tmp, strict=1, options=creation_opt)
        outDataset_cog = None
        ds_tmp = None
        tile_vec = None


def create_acq_id_image(outname, ref_tif, snap_datamasks, src_ids, extent,
                        epsg, driver, creation_opt, overviews, dst_nodata):
    """
    Creation of the Acquisition ID image.

    Parameters
    ----------
    outname: str
        Full path to the output data mask file.
    ref_tif: str
        Full path to any GeoTIFF file of the NRB product.
    snap_datamasks: list[str]
        List of raster datamask files covering the footprint of each source SLC scene that overlaps with the current
        MGRS tile.
    src_ids: list[pyroSAR.drivers.ID]
        List of :class:`~pyroSAR.drivers.ID` objects of all source SLC scenes that overlap with the current MGRS tile.
    extent: dict
        Spatial extent of the MGRS tile, derived from a :class:`~spatialist.vector.Vector` object.
    epsg: int
        The CRS used for the NRB product; provided as an EPSG code.
    driver: str
        GDAL driver to use for raster file creation.
    creation_opt: list[str]
        GDAL creation options to use for raster file creation. Should match specified GDAL driver.
    overviews: list[int]
        Internal overview levels to be created for each raster file.
    dst_nodata: int or str
        Nodata value to write to the output raster.
    """
    print(outname)
    src_scenes = [sid.scene for sid in src_ids]
    # If there are two source scenes, make sure that the order of acquisitions in all lists is correct!
    if len(src_scenes) > 1:
        if not len(src_scenes) == 2 and len(snap_datamasks) == 2:
            raise RuntimeError('expected lists `src_scenes` and `valid_mask_list` to be of length 2; length is '
                               '{} and {} respectively'.format(len(src_scenes), len(snap_datamasks)))
        starts_src = [datetime.strptime(identify(f).start, '%Y%m%dT%H%M%S') for f in src_scenes]
        start_valid = [datetime.strptime(re.search('[0-9]{8}T[0-9]{6}', os.path.basename(f)).group(), '%Y%m%dT%H%M%S')
                       for f in snap_datamasks]
        if starts_src[0] > starts_src[1]:
            src_scenes.reverse()
            starts_src.reverse()
        if start_valid[0] != starts_src[0]:
            snap_datamasks.reverse()
        if start_valid[0] != starts_src[0]:
            raise RuntimeError('failed to match order of lists `src_scenes` and `valid_mask_list`')
    
    tile_bounds = [extent['xmin'], extent['ymin'], extent['xmax'], extent['ymax']]
    
    arr_list = []
    for dm in snap_datamasks:
        vrt_snap_valid = '/vsimem/' + os.path.dirname(outname) + 'mosaic.vrt'
        gdalbuildvrt(dm, vrt_snap_valid, options={'outputBounds': tile_bounds}, void=False)
        with bbox(extent, crs=epsg) as tile_vec:
            with Raster(vrt_snap_valid)[tile_vec] as vrt_ras:
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
