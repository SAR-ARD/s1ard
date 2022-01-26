import os
import sys
import logging
from datetime import datetime
from lxml import etree
import binascii
from time import gmtime, strftime
import numpy as np
from osgeo import gdal
from spatialist import gdalbuildvrt, Raster
from spatialist.vector import bbox
from pyroSAR import identify, finder
from pyroSAR.ancillary import groupbyTime, seconds


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
    #consecutive_acq = [x for x in groupbyTime(images=selection, function=seconds, time=60) if isinstance(x, list)]
    
    #list_out = []
    #for scene in selection:
    #    sid = identify(scene)
    #    if len(finder(processdir, [sid.start], regex=True, recursive=False)) < 3 \
    #            and not any(scene in sublist for sublist in consecutive_acq):
    #        list_out.append(scene)
    
    #return list_out
    
    for scene in selection:
        list_all = finder(processdir, [identify(scene).start], regex=True, recursive=False)
        exclude = ['_NEBZ', '_NESZ', '_NEGZ']
        list_out = [scene for scene in selection if len([item for item in list_all
                                                         if not any(ex in item for ex in exclude)]) < 3]
    n_skipped = len(selection)-len(list_out)
    if n_skipped > 0:
        print("### {} scenes skipped, because they have already been processed.\n".format(n_skipped))
    print("### {} scenes in final selection for processing.".format(len(list_out)))
    return list_out


def modify_data_mask(dm_path, extent, epsg, driver, creation_opt, overviews, multilayer=False,
                     wbm=False, wbm_path=None):
    """
    Modifies the Data Mask file.
    
    Parameters
    ----------
    dm_path: str
        Path to the data mask.
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
        Should individual masks be written into seperate bands, creating a multi-level raster file? Default is False.
    wbm: bool, optional
        Include 'ocean water' information from an external Water Body Mask? Default is False.
    wbm_path: str, optional
        Path to the external Water Body Mask file.
    
    Returns
    -------
    None
    """
    ml_cog_out = dm_path.replace('.tif', '_2.tif')
    nodata = 255
    
    dm_bands = {1: {'arr_val': 0,
                    'name': 'not layover, nor shadow'},
                2: {'arr_val': 1,
                    'name': 'layover'},
                3: {'arr_val': 2,
                    'name': 'shadow'},
                4: {'arr_val': 4,
                    'name': 'ocean water'}}
    
    if wbm:
        assert os.path.isfile(wbm_path), "Water Body Mask '{}' not found".format(wbm_path)
        
        with bbox(extent, crs=epsg) as vec:
            with Raster(wbm_path)[vec] as ras_wbm:
                with Raster(dm_path) as ras_dm:
                    rows = ras_dm.rows
                    cols = ras_dm.cols
                    geotrans = ras_dm.raster.GetGeoTransform()
                    proj = ras_dm.raster.GetProjection()
                    
                    wbm_arr = ras_wbm.array()
                    dm_arr = ras_dm.array()
                    wbm_arr = np.where((wbm_arr == 1), 4, np.nan)
                    wbm_arr[np.isnan(wbm_arr)] = dm_arr[np.isnan(wbm_arr)]
                    
                    mask_arr = wbm_arr
                    del dm_arr
                    del wbm_arr
            vec = None
    else:
        with Raster(dm_path) as ras_dm:
            mask_arr = ras_dm.array()
            rows = ras_dm.rows
            cols = ras_dm.cols
            geotrans = ras_dm.raster.GetGeoTransform()
            proj = ras_dm.raster.GetProjection()
        dm_bands.pop(4)
    
    # MULTI-LAYER COG
    if multilayer:
        outname_tmp = '/vsimem/' + os.path.basename(ml_cog_out) + '.vrt'
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
            
            if arr_val == 0:
                arr = np.isnan(mask_arr)
            elif arr_val in [1, 2]:
                arr = np.where(((mask_arr == arr_val) | (mask_arr == 3)), 1, 0)
            elif arr_val == 4:
                arr = np.where(mask_arr == 4, 1, 0)
            
            arr[np.isnan(arr)] = 0
            arr = arr.astype('uint8')
            band.WriteArray(arr)
            del arr
            band.SetNoDataValue(nodata)
            band.SetDescription(b_name)
            band.FlushCache()
            band = None
        
        ds_tmp.SetMetadataItem('TIFFTAG_DATETIME', strftime('%Y:%m:%d %H:%M:%S', gmtime()))
        ds_tmp.BuildOverviews('AVERAGE', [2, 4, 8, 16, 32])
        outDataset_cog = gdal.GetDriverByName(driver).CreateCopy(ml_cog_out, ds_tmp,
                                                                 strict=1, options=creation_opt)
        outDataset_cog = None
        ds_tmp = None
    
    # SINGLE-LAYER COG
    with Raster(dm_path) as ras_dm:
        ras_dm.write(dm_path, format=driver, array=mask_arr.astype('uint8'), nodata=nodata, overwrite=True,
                     overviews=overviews, options=creation_opt)
    del mask_arr


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
