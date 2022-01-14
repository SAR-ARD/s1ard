import os
import sys
import logging
from datetime import datetime
from lxml import etree
import binascii
from spatialist import gdalbuildvrt
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
    
    list_out = [scene for scene in selection if len(finder(processdir, [identify(scene).start], regex=True,
                                                           recursive=False)) < 3]
    print("### {} scenes skipped, because they have already been processed.\n"
          "### {} scenes in final selection for processing.".format(len(selection)-len(list_out), len(list_out)))
    if len(list_out) == 0:
        sys.exit(0)
    else:
        return list_out


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
    ext_dem_file = {config['ext_dem_file']}
    ext_wbm_file = {config['ext_wbm_file']}
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
