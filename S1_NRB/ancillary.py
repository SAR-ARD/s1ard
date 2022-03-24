import os
import re
import sys
import logging
import binascii
from time import gmtime, strftime
from copy import deepcopy
from datetime import datetime, timezone
from lxml import etree
import numpy as np
from osgeo import gdal
from scipy.interpolate import griddata
import spatialist
from spatialist import gdalbuildvrt, Raster, bbox
import pyroSAR
from pyroSAR import identify, examine
import S1_NRB
from S1_NRB.metadata.extract import get_uid_sid, etree_from_sid, find_in_annotation


def create_vrt(src, dst, fun, relpaths=False, scale=None, offset=None, options=None,
               overviews=None, overview_resampling=None):
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
        A PixelFunctionType that should be applied on the fly when opening the VRT file. The function is applied to a
        band that derives its pixel information from the source bands. A list of possible options can be found here:
        https://gdal.org/drivers/raster/vrt.html#default-pixel-functions
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
    options: dict, optional
        Additional parameters passed to gdal.BuildVRT. For possible options see:
        https://gdal.org/python/osgeo.gdal-module.html#BuildVRTOptions
    overviews: list[int], optional
        Internal overview levels to be created for each raster file.
    overview_resampling: str, optional
        Resampling method for overview levels.
    
    Returns
    -------
    None
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
    Creates the RGB VRT file.
    
    Parameters
    ----------
    outname: str
        Full path to the output RGB VRT file.
    infiles: list[str]
        A list of paths.
    overviews: list[int]
        Internal overview levels to be defined for the created VRT file.
    overview_resampling: str
        Resampling method applied to overview pyramids.
    
    Returns
    -------
    None
    """
    print(outname)
    
    # make sure order is right and VV polarization is first
    paths_reorder = []
    for i, f in enumerate(infiles):
        pol = re.search('[hv]{2}', os.path.basename(f)).group()
        if i == 1 and pol == 'vv':
            paths_reorder.append(f)
            paths_reorder.append(infiles[0])
            infiles = paths_reorder
    
    # create VRT and change its content
    gdalbuildvrt(src=infiles, dst=outname, options={'separate': True})
    
    tree = etree.parse(outname)
    root = tree.getroot()
    bands = tree.findall('VRTRasterBand')
    new_band = etree.SubElement(root, 'VRTRasterBand',
                                attrib={'dataType': 'Float32', 'band': '3', 'subClass': 'VRTDerivedRasterBand'})
    vrt_nodata = etree.SubElement(new_band, 'NoDataValue')
    vrt_nodata.text = 'nan'
    new_band.insert(1, deepcopy(bands[0].find('ComplexSource')))
    new_band.insert(2, deepcopy(bands[1].find('ComplexSource')))
    # pxfun_type = etree.SubElement(new_band, 'PixelFunctionType')
    # pxfun_type.text = 'diff'
    pxfun_language = etree.SubElement(new_band, 'PixelFunctionLanguage')
    pxfun_language.text = 'Python'
    pxfun_type = etree.SubElement(new_band, 'PixelFunctionType')
    pxfun_type.text = 'div'
    pxfun_code = etree.SubElement(new_band, 'PixelFunctionCode')
    pxfun_code.text = etree.CDATA("""
import numpy as np
def div(in_ar, out_ar, xoff, yoff, xsize, ysize, raster_xsize, raster_ysize, buf_radius, gt, **kwargs):
    np.divide(in_ar[0], in_ar[1], out=out_ar, where=in_ar[1]!=0, dtype='float32')
    """)
    
    bands = tree.findall('VRTRasterBand')
    for band, col in zip(bands, ['Red', 'Green', 'Blue']):
        color = etree.Element('ColorInterp')
        color.text = col
        band.insert(1, color)
        if band.attrib['band'] in ['1', '2']:
            ndv = band.find('NoDataValue')
            band.remove(ndv)
    
    ovr = etree.SubElement(root, 'OverviewList', attrib={'resampling': overview_resampling.lower()})
    ov = str(overviews)
    for x in ['[', ']', ',']:
        ov = ov.replace(x, '')
    ovr.text = ov
    
    etree.indent(root)
    tree.write(outname, pretty_print=True, xml_declaration=False, encoding='utf-8')


def generate_unique_id(encoded_str):
    """
    Returns a unique product identifier as a hexa-decimal string generated from the time of execution in isoformat.
    The CRC-16 algorithm used to compute the unique identifier is CRC-CCITT (0xFFFF).
    
    Parameters
    ----------
    encoded_str: bytes
        A string that should be used to generate a unique id from. The string needs to be encoded; e.g.:
        `'abc'.encode()`
    
    Returns
    -------
    p_id: str
        The unique product identifier.
    """
    crc = binascii.crc_hqx(encoded_str, 0xffff)
    p_id = '{:04X}'.format(crc & 0xffff)
    
    return p_id


def calc_product_start_stop(src_scenes, extent, epsg):
    """
    Calculates the start and stop times of the NRB product.
    The geolocation grid points including their azimuth time information are extracted first from the metadata of each
    source SLC. These grid points are then used to interpolate the azimuth time for the lower right and upper left
    (ascending) or upper right and lower left (descending) corners of the MGRS tile extent.
    
    Parameters
    ----------
    src_scenes: list[str]
        A list of paths pointing to the source scenes of the NRB product.
    extent: dict
        Spatial extent of the MGRS tile, derived from a `spatialist.vector.Vector` object.
    epsg: int
        The CRS used for the NRB product; provided as an EPSG code.
    
    Returns
    -------
    start: str
        Start time of the NRB product formatted as %Y%m%dT%H%M%S in UTC.
    stop: str
        Stop time of the NRB product formatted as %Y%m%dT%H%M%S in UTC.
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
    for i in range(len(src_scenes)):
        uid, sid = get_uid_sid(src_scenes[i])
        slc_dict[uid] = etree_from_sid(sid)
        slc_dict[uid]['sid'] = sid
    
    uids = list(slc_dict.keys())
    swaths = list(slc_dict[uids[0]]['annotation'].keys())
    
    for uid in uids:
        t = find_in_annotation(annotation_dict=slc_dict[uid]['annotation'], pattern='.//geolocationGridPoint/azimuthTime')
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


def create_data_mask(outname, valid_mask_list, snap_files, extent, epsg, driver, creation_opt, overviews,
                     overview_resampling, wbm=None):
    """
    Creates the Data Mask file.
    
    Parameters
    ----------
    outname: str
        Full path to the output data mask file.
    valid_mask_list: list[str]
        A list of paths pointing to the datamask_ras files that intersect with the current MGRS tile.
    snap_files: list[str]
        A list of paths pointing to the SNAP processed datasets of the NRB product.
    extent: dict
        Spatial extent of the MGRS tile, derived from a `spatialist.vector.Vector` object.
    epsg: int
        The CRS used for the NRB product; provided as an EPSG code.
    driver: str
        GDAL driver to use for raster file creation.
    creation_opt: list[str]
        GDAL creation options to use for raster file creation. Should match specified GDAL driver.
    overviews: list[int]
        Internal overview levels to be created for each raster file.
    overview_resampling: str
        Resampling method for overview levels.
    wbm: str, optional
        Path to a water body mask file with the dimensions of an MGRS tile.
    
    Returns
    -------
    None
    """
    print(outname)
    out_nodata = 255
    
    pols = [pol for pol in set([re.search('[VH]{2}', os.path.basename(x)).group() for x in snap_files if
                                re.search('[VH]{2}', os.path.basename(x)) is not None])]
    pattern = pols[0] + '_gamma0-rtc'
    snap_gamma0 = [x for x in snap_files if re.search(pattern, os.path.basename(x))]
    snap_ls_mask = [x for x in snap_files if re.search('layoverShadowMask', os.path.basename(x))]
    
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
    gdalbuildvrt(valid_mask_list, vrt_snap_valid, options={'outputBounds': tile_bounds}, void=False)
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
                    out_arr[np.isnan(arr_snap_valid)] = out_nodata
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
            arr[out_arr == out_nodata] = out_nodata
            if arr_val == 0:
                arr[out_arr == 0] = 1
            elif arr_val in [1, 2]:
                arr[(out_arr == arr_val) | (out_arr == 3)] = 1
            elif arr_val == 4:
                arr[out_arr == 4] = 1
            
            arr = arr.astype('uint8')
            band.WriteArray(arr)
            band.SetNoDataValue(out_nodata)
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


def create_acq_id_image(ref_tif, valid_mask_list, src_scenes, extent, epsg, driver, creation_opt, overviews):
    """
    Creation of the acquisition ID image described in CARD4L NRB 2.8
    
    Parameters
    ----------
    ref_tif: str
        Full path to any GeoTIFF file of the NRB product.
    valid_mask_list: list[str]
        A list of paths pointing to the datamask_ras files that intersect with the current MGRS tile.
    src_scenes: list[str]
        A list of paths pointing to the source scenes of the NRB product.
    extent: dict
        Spatial extent of the MGRS tile, derived from a `spatialist.vector.Vector` object.
    epsg: int
        The CRS used for the NRB product; provided as an EPSG code.
    driver: str
        GDAL driver to use for raster file creation.
    creation_opt: list[str]
        GDAL creation options to use for raster file creation. Should match specified GDAL driver.
    overviews: list[int]
        Internal overview levels to be created for each raster file.
    
    Returns
    -------
    None
    """
    outname = ref_tif.replace('-gs.tif', '-id.tif')
    print(outname)
    out_nodata = 255
    
    # If there are two source scenes, make sure that the order of acquisitions in all lists is correct!
    if len(src_scenes) > 1:
        if not len(src_scenes) == 2 and len(valid_mask_list) == 2:
            raise RuntimeError('expected lists `src_scenes` and `valid_mask_list` to be of length 2; length is '
                               '{} and {} respectively'.format(len(src_scenes), len(valid_mask_list)))
        starts_src = [datetime.strptime(identify(f).start, '%Y%m%dT%H%M%S') for f in src_scenes]
        start_valid = [datetime.strptime(re.search('[0-9]{8}T[0-9]{6}', os.path.basename(f)).group(), '%Y%m%dT%H%M%S')
                       for f in valid_mask_list]
        if starts_src[0] > starts_src[1]:
            src_scenes.reverse()
            starts_src.reverse()
        if start_valid[0] != starts_src[0]:
            valid_mask_list.reverse()
        if start_valid[0] != starts_src[0]:
            raise RuntimeError('failed to match order of lists `src_scenes` and `valid_mask_list`')
    
    tile_bounds = [extent['xmin'], extent['ymin'], extent['xmax'], extent['ymax']]
    
    arr_list = []
    for file in valid_mask_list:
        vrt_snap_valid = '/vsimem/' + os.path.dirname(outname) + 'mosaic.vrt'
        gdalbuildvrt(file, vrt_snap_valid, options={'outputBounds': tile_bounds}, void=False)
        with bbox(extent, crs=epsg) as tile_vec:
            with Raster(vrt_snap_valid)[tile_vec] as vrt_ras:
                vrt_arr = vrt_ras.array()
                arr_list.append(vrt_arr)
                del vrt_arr
            tile_vec = None
    
    src_scenes_clean = [os.path.basename(src).replace('.zip', '').replace('.SAFE', '') for src in src_scenes]
    tag = '{{"{src1}": 1}}'.format(src1=src_scenes_clean[0])
    out_arr = np.full(arr_list[0].shape, out_nodata)
    out_arr[arr_list[0] == 1] = 1
    if len(arr_list) == 2:
        out_arr[arr_list[1] == 1] = 2
        tag = '{{"{src1}": 1, "{src2}": 2}}'.format(src1=src_scenes_clean[0], src2=src_scenes_clean[1])
    
    creation_opt.append('TIFFTAG_IMAGEDESCRIPTION={}'.format(tag))
    with Raster(ref_tif) as ref_ras:
        ref_ras.write(outname, format=driver, array=out_arr.astype('uint8'), nodata=out_nodata, overwrite=True,
                      overviews=overviews, options=creation_opt)


def get_max_ext(geometries, buffer=None):
    """
    Gets the maximum extent from a list of geometries.
    
    Parameters
    ----------
    geometries: list[spatialist.vector.Vector objects]
        List of vector geometries.
    buffer: float, optional
        The buffer in degrees to add to the extent.
    Returns
    -------
    max_ext: dict
        The maximum extent of the selected vector geometries including the chosen buffer.
    """
    max_ext = {}
    for geo in geometries:
        if len(max_ext.keys()) == 0:
            max_ext = geo.extent
        else:
            for key in ['xmin', 'ymin']:
                if geo.extent[key] < max_ext[key]:
                    max_ext[key] = geo.extent[key]
            for key in ['xmax', 'ymax']:
                if geo.extent[key] > max_ext[key]:
                    max_ext[key] = geo.extent[key]
    max_ext = dict(max_ext)
    if buffer is not None:
        max_ext['xmin'] -= buffer
        max_ext['xmax'] += buffer
        max_ext['ymin'] -= buffer
        max_ext['ymax'] += buffer
    return max_ext


def set_logging(config, debug=False):
    """
    Set logging for the current process.
    
    Parameters
    ----------
    config: dict
        Dictionary of the parsed config parameters for the current process.
    debug: bool, optional
        Set pyroSAR logging level to DEBUG? Default is False.
    
    Returns
    -------
    log_local: logging.Logger
        The log handler for the current process.
    """
    # pyroSAR logging as sys.stdout
    log_pyro = logging.getLogger('pyroSAR')
    if debug:
        log_pyro.setLevel(logging.DEBUG)
    else:
        log_pyro.setLevel(logging.INFO)
    sh = logging.StreamHandler(sys.stdout)
    log_pyro.addHandler(sh)
    
    # NRB logging in logfile
    now = datetime.now().strftime('%Y%m%dT%H%M')
    log_local = logging.getLogger(__name__)
    log_local.setLevel(logging.DEBUG)
    log_file = os.path.join(config['log_dir'], f"{now}_process.log")
    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    fh = logging.FileHandler(filename=log_file, mode='a')
    log_local.addHandler(fh)
    
    # Add header first with simple formatting
    form_simple = logging.Formatter("%(message)s")
    fh.setFormatter(form_simple)
    _log_process_config(logger=log_local, config=config)
    
    # Use normal formatting from here on out
    form = logging.Formatter("[%(asctime)s] [%(levelname)8s] %(message)s")
    fh.setFormatter(form)
    
    return log_local


def _log_process_config(logger, config):
    """
    Adds a header to the logfile, which includes information about the current processing configuration.
    
    Parameters
    ----------
    logger: Logger
        The logger to which the header is added to.
    config: dict
        Dictionary of the parsed config parameters for the current process.
    
    Returns
    -------
    None
    """
    try:
        core = examine.ExamineSnap().get_version('core')
        s1tbx = examine.ExamineSnap().get_version('s1tbx')
        snap_core = f"{core['version']} | {core['date']}"
        snap_s1tbx = f"{s1tbx['version']} | {s1tbx['date']}"
    except RuntimeError:
        snap_core = 'unknown'
        snap_s1tbx = 'unknown'
    
    header = f"""
    ====================================================================================================================
    PROCESSING CONFIGURATION
    
    mode = {config['mode']}
    aoi_tiles = {config['aoi_tiles']}
    aoi_geometry = {config['aoi_geometry']}
    mindate = {config['mindate'].isoformat()}
    maxdate = {config['maxdate'].isoformat()}
    acq_mode = {config['acq_mode']}
    
    work_dir = {config['work_dir']}
    scene_dir = {config['scene_dir']}
    rtc_dir = {config['rtc_dir']}
    tmp_dir = {config['tmp_dir']}
    nrb_dir = {config['nrb_dir']}
    dem_dir = {config['dem_dir']}
    wbm_dir = {config['wbm_dir']}
    log_dir = {config['log_dir']}
    db_file = {config['db_file']}
    kml_file = {config['kml_file']}
    dem_type = {config.get('dem_type')}
    
    ====================================================================================================================
    SOFTWARE
    
    S1_NRB: {S1_NRB.__version__}
    snap-core: {snap_core}
    snap-s1tbx: {snap_s1tbx}
    python: {sys.version}
    python-pyroSAR: {pyroSAR.__version__}
    python-spatialist: {spatialist.__version__}
    python-GDAL: {gdal.__version__}
    gdal_threads = {config.get('gdal_threads')}
    
    ====================================================================================================================
    """
    logger.info(header)
