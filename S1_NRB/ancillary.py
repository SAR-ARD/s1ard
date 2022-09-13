import os
import sys
import logging
import binascii
from datetime import datetime
from osgeo import gdal
import spatialist
import pyroSAR
from pyroSAR import examine
import S1_NRB


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


def get_max_ext(geometries, buffer=None):
    """
    Gets the maximum extent from a list of geometries.
    
    Parameters
    ----------
    geometries: list[spatialist.vector.Vector]
        List of :class:`~spatialist.vector.Vector` geometries.
    buffer: float, optional
        The buffer in degrees to add to the extent.
    Returns
    -------
    max_ext: dict
        The maximum extent of the selected :class:`~spatialist.vector.Vector` geometries including the chosen buffer.
    """
    max_ext = {}
    for geo in geometries:
        if len(max_ext.keys()) == 0:
            max_ext = geo.extent
        else:
            ext = geo.extent
            for key in ['xmin', 'ymin']:
                if ext[key] < max_ext[key]:
                    max_ext[key] = ext[key]
            for key in ['xmax', 'ymax']:
                if ext[key] > max_ext[key]:
                    max_ext[key] = ext[key]
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
    gdal_threads = {config.get('gdal_threads')}
    
    ====================================================================================================================
    SOFTWARE
    
    S1_NRB: {S1_NRB.__version__}
    snap-core: {snap_core}
    snap-s1tbx: {snap_s1tbx}
    python: {sys.version}
    python-pyroSAR: {pyroSAR.__version__}
    python-spatialist: {spatialist.__version__}
    python-GDAL: {gdal.__version__}
    
    ====================================================================================================================
    """
    logger.info(header)


def log(handler, mode, proc_step, scenes, epsg, msg):
    """
    Format and handle log messages during processing.
    
    Parameters
    ----------
    handler: logging.Logger
        The log handler for the current process.
    mode: str
        One of ['info', 'warning', 'exception']. Calls the respective logging helper function. E.g., `handler.info()`.
    proc_step: str
        The processing step for which the message is logged.
    scenes: str or list[str]
        Scenes that are currently being processed.
    epsg: int
        The coordinate reference system as an EPSG code.
    msg: Any
        The massage that should be logged.
    """
    proc_step = proc_step.zfill(7).replace('0', ' ')
    message = '[{proc_step}] -- {scenes} [{epsg}] -- {msg}'
    message = message.format(proc_step=proc_step, scenes=scenes, epsg=epsg, msg=msg)
    if mode == 'info':
        handler.info(message)
    elif mode == 'warning':
        handler.warning(message)
    elif mode == 'exception':
        handler.exception(message)
    else:
        raise RuntimeError('log mode {} is not supported'.format(mode))
