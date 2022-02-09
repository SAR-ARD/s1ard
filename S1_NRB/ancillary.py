import os
import re
import sys
import logging
from datetime import datetime
from lxml import etree
import binascii
from time import gmtime, strftime
import numpy as np
from osgeo import gdal
from scipy.interpolate import griddata
from spatialist import gdalbuildvrt, Raster, bbox
from pyroSAR import identify, finder
from pyroSAR.ancillary import groupbyTime, seconds

from S1_NRB.metadata.extract import get_uid_sid, etree_from_sid, find_in_annotation


def vrt_pixfun(src, dst, fun, scale=None, offset=None, options=None):
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
    scale: int, optional
         The scale that should be applied when computing “real” pixel values from scaled pixel values on a raster band.
    offset: float, optional
        The offset that should be applied when computing “real” pixel values from scaled pixel values on a raster band.
    options: dict, optional
        Additional parameters passed to gdal.BuildVRT. For possible options see:
        https://gdal.org/python/osgeo.gdal-module.html#BuildVRTOptions
    
    Returns
    -------
    None
    """
    gdalbuildvrt(src=src, dst=dst, options=options)
    tree = etree.parse(dst)
    band = tree.find('VRTRasterBand')
    band.attrib['subClass'] = 'VRTDerivedRasterBand'
    pixfun = etree.SubElement(band, 'PixelFunctionType')
    pixfun.text = fun
    if scale is not None:
        sc = etree.SubElement(band, 'Scale')
        sc.text = str(scale)
    if offset is not None:
        off = etree.SubElement(band, 'Offset')
        off.text = str(offset)
    complexSrc = band.find('ComplexSource')
    nodata = complexSrc.find('NODATA')
    nodata.text = '0'
    tree.write(dst, pretty_print=True, xml_declaration=False, encoding='utf-8')


def vrt_relpath(vrt):
    """
    Converts the paths to annotation datasets in a VRT file to relative paths.
    
    Parameters
    ----------
    vrt: str
        Path to the VRT file that should be updated.
    
    Returns
    -------
    None
    """
    tree = etree.parse(vrt)
    test = tree.xpath('//SourceFilename[@relativeToVRT="0"]')[0]
    repl = '../annotation/' + os.path.basename(test.text.replace('\\', '\\\\'))
    test.text = repl
    test.attrib['relativeToVRT'] = '1'
    tree.write(vrt, pretty_print=True, xml_declaration=False, encoding='utf-8')


def generate_product_id():
    """
    Returns a unique product identifier as a hexa-decimal string generated from the time of execution in isoformat.
    The CRC-16 algorithm used to compute the unique identifier is CRC-CCITT (0xFFFF).
    
    Returns
    -------
    p_id: str
        The unique product identifier.
    t: datetime.datetime
        The datetime object used to generate the unique product identifier from.
    """
    t = datetime.now()
    tie = t.isoformat().encode()
    crc = binascii.crc_hqx(tie, 0xffff)
    p_id = '{:04X}'.format(crc & 0xffff)
    
    return p_id, t


def filter_selection(selection, processdir):
    """
    Filters a selection of scenes that should be processed. Returns only those scenes that have not been (successfully)
    processed before. If more than 3 files belonging to a scene are found, it is assumed that it was successfully
    processed before and filtered from the selection.
    
    Parameters
    ----------
    selection: list[str]
        The selection of scenes that should be filtered.
    processdir: str
        The directory to scan for already processed scenes.
    
    Returns
    -------
    list_out: list[str]
        The filtered selection of scenes.
    """
    list_out = []
    for scene in selection:
        list_processed = finder(processdir, [identify(scene).start], regex=True, recursive=False)
        exclude = ['_NEBZ', '_NESZ', '_NEGZ']
        if len([item for item in list_processed if not any(ex in item for ex in exclude)]) < 3:
            list_out.append(scene)
    
    n_skipped = len(selection)-len(list_out)
    if n_skipped > 0:
        print("### {} scenes skipped, because they have already been processed.\n".format(n_skipped))
    print("### {} scenes in final selection for processing.".format(len(list_out)))
    return list_out


def calc_product_start_stop(src_scenes, extent, epsg):
    """
    Calculates the start and stop times of the current product. The geolocation grid points including their azimuth time
    information are extracted first from the metadata of each source SLC. These grid points are then used to interpolate
    the azimuth time for the lower right and upper left (Ascending) or upper right and lower left (Descending) corners
    of the MGRS tile of the current product.
    
    Parameters
    ----------
    src_scenes: list[str]
        A list of paths pointing to the source scenes of the product.
    extent: dict
        Spatial extent of the MGRS tile, derived from a `spatialist.vector.Vector` object.
    epsg: int
        The CRS used for the NRB product; provided as an EPSG code.
    
    Returns
    -------
    start: str
        Start time of the current product formatted as %Y%m%dT%H%M%S.
    stop: str
        Stop time of the current product formatted as %Y%m%dT%H%M%S.
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


def create_data_mask(outname, valid_mask_list, src_files, extent, epsg, driver, creation_opt, overviews,
                     multilayer=False, wbm=False, wbm_path=None):
    """
    Creates the Data Mask file.
    
    Parameters
    ----------
    outname: str
        Full path to the output data mask file.
    valid_mask_list: list[str]
        A list of paths pointing to the datamask_ras files that intersect with the current MGRS tile.
    src_files: list[str]
        A list of paths pointing to the SNAP processed datasets of the product.
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
    multilayer: bool, optional
        Should individual masks be written into separate bands, creating a multi-level raster file? Default is False.
    wbm: bool, optional
        Include 'ocean water' information from an external Water Body Mask? Default is False.
    wbm_path: str, optional
        Path to the external Water Body Mask file. Ignored if `wbm=False`.
    
    Returns
    -------
    None
    """
    print(outname)
    outname_ml = outname.replace('.tif', '_2.tif')
    out_nodata = 255
    
    pols = [pol for pol in set([re.search('[VH]{2}', os.path.basename(x)).group() for x in src_files if
                                re.search('[VH]{2}', os.path.basename(x)) is not None])]
    pattern = pols[0] + '_gamma0-rtc'
    snap_gamma0 = [x for x in src_files if re.search(pattern, os.path.basename(x))]
    snap_ls_mask = [x for x in src_files if re.search('layoverShadowMask', os.path.basename(x))]
    
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
            
            # Add Water Body Mask if wbm=True
            if wbm:
                with Raster(wbm_path)[tile_vec] as ras_wbm:
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
        
        # SINGLE-LAYER COG
        ras_snap_ls.write(outname, format=driver,
                          array=out_arr.astype('uint8'), nodata=out_nodata, overwrite=True, overviews=overviews,
                          options=creation_opt)
        
        # MULTI-LAYER COG
        if multilayer:
            outname_tmp = '/vsimem/' + os.path.basename(outname_ml) + '.vrt'
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
                
                arr = np.full((rows, cols), out_nodata)
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
            ds_tmp.BuildOverviews('AVERAGE', [2, 4, 8, 16, 32])
            outDataset_cog = gdal.GetDriverByName(driver).CreateCopy(outname_ml, ds_tmp,
                                                                     strict=1, options=creation_opt)
            outDataset_cog = None
            ds_tmp = None
        tile_vec = None


def create_acq_id_image(ref_tif, valid_mask_list, src_scenes, extent, epsg, driver, creation_opt, overviews):
    """
    Creation of the acquisition ID image described in CARD4L NRB 2.8
    
    Parameters
    ----------
    ref_tif: str
        Full path to a reference GeoTIFF file of the product.
    valid_mask_list: list[str]
        A list of paths pointing to the datamask_ras files that intersect with the current MGRS tile.
    src_scenes: list[str]
        A list of paths pointing to the source scenes of the product.
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
    
    # If there are two source scenes, make sure that the order in the relevant lists is correct!
    if len(src_scenes) == 2:
        starts = [datetime.strptime(identify(f).start, '%Y%m%dT%H%M%S') for f in src_scenes]
        if starts[0] > starts[1]:
            src_scenes_new = [src_scenes[1]]
            src_scenes_new.append(src_scenes[0])
            src_scenes = src_scenes_new
            starts = [identify(f).start for f in src_scenes]
        start_valid = [datetime.strptime(re.search('[0-9]{8}T[0-9]{6}', os.path.basename(f)).group(),
                                         '%Y%m%dT%H%M%S') for f in valid_mask_list]
        if start_valid[0] != starts[0]:
            valid_mask_list_new = [valid_mask_list[1]]
            valid_mask_list_new.append(valid_mask_list[0])
            valid_mask_list = valid_mask_list_new
    
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
    tag = '\n1: {src1}'.format(src1=src_scenes_clean[0])
    out_arr = np.full(arr_list[0].shape, out_nodata)
    if len(arr_list) == 2:
        out_arr[arr_list[1] == 1] = 2
        tag = '\n1: {src1}\n2: {src2}'.format(src1=src_scenes_clean[0], src2=src_scenes_clean[1])
    
    out_arr[arr_list[0] == 1] = 1
    creation_opt.append('TIFFTAG_IMAGEDESCRIPTION={}'.format(tag))
    
    with Raster(ref_tif) as ref_ras:
        ref_ras.write(outname, format=driver, array=out_arr.astype('uint8'), nodata=out_nodata, overwrite=True,
                      overviews=overviews, options=creation_opt)


def set_logging(config):
    """
    Set logging for the current process.
    
    Parameters
    ----------
    config: dict
        Dictionary of the parsed config parameters for the current process.
    
    Returns
    -------
    log_local: logging.Logger
        The log handler for the current process.
    """
    
    # pyroSAR logging as sys.stdout
    log_pyro = logging.getLogger('pyroSAR')
    log_pyro.setLevel(logging.INFO)
    sh = logging.StreamHandler(sys.stdout)
    log_pyro.addHandler(sh)
    
    # NRB logging in logfile
    now = datetime.now().strftime('%Y%m%dT%H%M')
    log_local = logging.getLogger(__name__)
    log_local.setLevel(logging.DEBUG)
    log_file = os.path.join(config['work_dir'], 'log', f"{now}_process.log")
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
    
    import pyroSAR
    from pyroSAR.examine import ExamineSnap
    from osgeo import gdal
    
    core = ExamineSnap().get_version('core')
    s1tbx = ExamineSnap().get_version('s1tbx')
    
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
    out_dir = {config['out_dir']}
    tmp_dir = {config['tmp_dir']}
    db_file = {config['db_file']}
    kml_file = {config['kml_file']}
    ext_dem_file = {config.get('ext_dem_file')}
    ext_wbm_file = {config.get('ext_wbm_file')}
    ====================================================================================================================
    SOFTWARE
    
    snap-core: {core['version']} | {core['date']}
    snap-s1tbx: {s1tbx['version']} | {s1tbx['date']}
    python: {sys.version}
    python-pyroSAR: {pyroSAR.__version__}
    python-GDAL: {gdal.__version__}
    
    ====================================================================================================================
    """
    
    logger.info(header)
