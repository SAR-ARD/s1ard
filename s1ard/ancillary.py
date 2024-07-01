import os
import sys
import logging
import requests
import binascii
from lxml import etree
from textwrap import dedent
from datetime import datetime, timedelta
from osgeo import gdal
import spatialist
from spatialist.vector import bbox, intersect
import pyroSAR
from pyroSAR.ancillary import Lock
from pyroSAR import examine, identify_many
import s1ard

log = logging.getLogger('s1ard')


def check_scene_consistency(scenes):
    """
    Check the consistency of a scene selection.
    The following pyroSAR object attributes must be the same:
    
     - sensor
     - acquisition_mode
     - product
     - frameNumber (data take ID)
    
    Parameters
    ----------
    scenes: list[str or pyroSAR.drivers.ID]

    Returns
    -------
    
    Raises
    ------
    RuntimeError
    """
    scenes = identify_many(scenes)
    for attr in ['sensor', 'acquisition_mode', 'product', 'frameNumber']:
        values = set([getattr(x, attr) for x in scenes])
        if not len(values) == 1:
            msg = f"scene selection differs in attribute '{attr}': {values}"
            raise RuntimeError(msg)


def check_spacing(spacing):
    """
    Check whether the spacing fits into the MGRS tile boundaries.

    Parameters
    ----------
    spacing: int or float
        the target pixel spacing in meters

    Returns
    -------

    """
    if 109800 % spacing != 0:
        raise RuntimeError(f'target spacing of {spacing} m does not align with MGRS tile size of 109800 m.')


def generate_unique_id(encoded_str):
    """
    
    Returns a unique product identifier as a hexadecimal string.
    The CRC-16 algorithm used to compute the unique identifier is CRC-CCITT (0xFFFF).
    
    Parameters
    ----------
    encoded_str: bytes
        A string that should be used to generate a unique id from. The string needs to be encoded; e.g.:
        ``'abc'.encode()``
    
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
    buffer: float or None
        The buffer in units of the geometries' CRS to add to the extent.
    
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
    debug: bool
        Set logging level to DEBUG?
    
    Returns
    -------
    logging.Logger
        The log handler for the current process.
    """
    level = logging.DEBUG if debug else logging.INFO
    
    logger = logging.getLogger('s1ard')
    logger.setLevel(level)
    
    log_format = "[%(asctime)s] [%(levelname)5s] %(message)s"
    formatter = logging.Formatter(fmt=log_format,
                                  datefmt='%Y-%m-%d %H:%M:%S')
    
    if config['logfile'] is not None:
        os.makedirs(os.path.dirname(config['logfile']), exist_ok=True)
        handler = logging.FileHandler(filename=config['logfile'], mode='a')
    else:
        handler = logging.StreamHandler(sys.stdout)
    logger.addHandler(handler)
    
    # Add header first with simple formatting
    formatter_simple = logging.Formatter("%(message)s")
    handler.setFormatter(formatter_simple)
    _log_process_config(logger=logger, config=config)
    
    # Use normal formatting from here on
    handler.setFormatter(formatter)
    
    # add pyroSAR logger
    log_pyro = logging.getLogger('pyroSAR')
    log_pyro.setLevel(level)
    log_pyro.addHandler(handler)
    
    return logger


def group_by_time(scenes, time=3):
    """
    Group scenes by their acquisition time difference.

    Parameters
    ----------
    scenes:list[pyroSAR.drivers.ID or str]
        a list of image names
    time: int or float
        a time difference in seconds by which to group the scenes.
        The default of 3 seconds incorporates the overlap between SLCs.

    Returns
    -------
    list[list[pyroSAR.drivers.ID]]
        a list of sub-lists containing the file names of the grouped scenes
    """
    # sort images by time stamp
    scenes = identify_many(scenes, sortkey='start')
    
    if len(scenes) < 2:
        return [scenes]
    
    groups = [[scenes[0]]]
    group = groups[0]
    
    for i in range(1, len(scenes)):
        start = datetime.strptime(scenes[i].start, '%Y%m%dT%H%M%S')
        stop_pred = datetime.strptime(scenes[i - 1].stop, '%Y%m%dT%H%M%S')
        diff = abs((stop_pred - start).total_seconds())
        if diff <= time:
            group.append(scenes[i])
        else:
            groups.append([scenes[i]])
            group = groups[-1]
    return groups


def _log_process_config(logger, config):
    """
    Adds a header to the logfile, which includes information about the current processing configuration.
    
    Parameters
    ----------
    logger: logging.Logger
        The logger to which the header is added to.
    config: dict
        Dictionary of the parsed config parameters for the current process.
    """
    try:
        snap_config = examine.ExamineSnap()
        core = snap_config.get_version('core')
        microwavetbx = snap_config.get_version('microwavetbx')
        snap_core = f"{core['version']} | {core['date']}"
        snap_microwavetbx = f"{microwavetbx['version']} | {microwavetbx['date']}"
    except RuntimeError:
        snap_core = 'unknown'
        snap_microwavetbx = 'unknown'
    
    header = f"""
    ====================================================================================================================
    PROCESSING CONFIGURATION
    
    mode                {config['mode']}
    aoi_tiles           {config['aoi_tiles']}
    aoi_geometry        {config['aoi_geometry']}
    scene               {config['scene']}
    mindate             {config['mindate'].isoformat()}
    maxdate             {config['maxdate'].isoformat()}
    date_strict         {config['date_strict']}
    sensor              {config['sensor']}
    acq_mode            {config['acq_mode']}
    product             {config['product']}
    datatake            {config['datatake']}
    measurement         {config.get('measurement')}
    annotation          {config.get('annotation')}
    dem_type            {config.get('dem_type')}
    etad                {config.get('etad')}
    
    work_dir            {config['work_dir']}
    sar_dir             {config['sar_dir']}
    tmp_dir             {config['tmp_dir']}
    ard_dir             {config['ard_dir']}
    wbm_dir             {config['wbm_dir']}
    etad_dir            {config['etad_dir']}
    scene_dir           {config['scene_dir']}
    logfile             {config['logfile']}
    db_file             {config['db_file']}
    stac_catalog        {config['stac_catalog']}
    stac_collections    {config['stac_collections']}
    gdal_threads        {config.get('gdal_threads')}
    snap_gpt_args       {config['snap_gpt_args']}
    
    ====================================================================================================================
    SOFTWARE
    
    s1ard               {s1ard.__version__}
    snap-core           {snap_core}
    snap-microwavetbx   {snap_microwavetbx}
    python              {sys.version}
    python-pyroSAR      {pyroSAR.__version__}
    python-spatialist   {spatialist.__version__}
    python-GDAL         {gdal.__version__}
    
    ====================================================================================================================
    """
    logger.info(dedent(header))


def vrt_add_overviews(vrt, overviews, resampling='AVERAGE'):
    """
    Add overviews to an existing VRT file.
    Existing overviews will be overwritten.

    Parameters
    ----------
    vrt: str
        the VRT file
    overviews: list[int]
         the overview levels
    resampling: str
        the overview resampling method

    Returns
    -------

    """
    tree = etree.parse(vrt)
    root = tree.getroot()
    ovr = root.find('OverviewList')
    if ovr is None:
        ovr = etree.SubElement(root, 'OverviewList')
    ovr.text = ' '.join([str(x) for x in overviews])
    ovr.attrib['resampling'] = resampling.lower()
    etree.indent(root)
    tree.write(vrt, pretty_print=True, xml_declaration=False, encoding='utf-8')


def buffer_min_overlap(geom1, geom2, percent=1):
    """
    Buffer a geometry to a minimum overlap with a second geometry.
    The geometry is iteratively buffered until the minimum overlap is reached.
    If the overlap of the input geometries is already larger than the defined
    threshold, a copy of the original geometry is returned.

    Parameters
    ----------
    geom1: spatialist.vector.Vector
        the geometry to be buffered
    geom2: spatialist.vector.Vector
        the reference geometry to intersect with
    percent: int or float
        the minimum overlap in percent of `geom1`

    Returns
    -------

    """
    geom2_area = geom2.getArea()
    ext = geom1.extent
    ext2 = ext.copy()
    xdist = ext['xmax'] - ext['xmin']
    ydist = ext['ymax'] - ext['ymin']
    buffer = 0
    overlap = 0
    while overlap < percent:
        xbuf = xdist * buffer / 100 / 2
        ybuf = ydist * buffer / 100 / 2
        ext2['xmin'] = ext['xmin'] - xbuf
        ext2['xmax'] = ext['xmax'] + xbuf
        ext2['ymin'] = ext['ymin'] - ybuf
        ext2['ymax'] = ext['ymax'] + ybuf
        with bbox(ext2, 4326) as geom3:
            ext3 = geom3.extent
            with intersect(geom2, geom3) as inter:
                inter_area = inter.getArea()
                overlap = inter_area / geom2_area * 100
        buffer += 1
    return bbox(ext3, 4326)


def buffer_time(start, stop, **kwargs):
    """
    Time range buffering
    
    Parameters
    ----------
    start: str
        the start time in format '%Y%m%dT%H%M%S'
    stop: str
        the stop time in format '%Y%m%dT%H%M%S'
    kwargs
        time arguments passed to :func:`datetime.timedelta`

    Returns
    -------

    """
    f = '%Y%m%dT%H%M%S'
    td = timedelta(**kwargs)
    start = datetime.strptime(start, f) - td
    start = datetime.strftime(start, f)
    stop = datetime.strptime(stop, f) + td
    stop = datetime.strftime(stop, f)
    return start, stop


def get_kml():
    """
    Download the S2 grid KML file to ~/s1ard and return the file path.
    
    Returns
    -------
    str
        the path to the KML file
    """
    remote = ('https://sentinel.esa.int/documents/247904/1955685/'
              'S2A_OPER_GIP_TILPAR_MPC__20151209T095117_V20150622T000000_21000101T000000_B00.kml')
    local_path = os.path.join(os.path.expanduser('~'), 's1ard')
    os.makedirs(local_path, exist_ok=True)
    local = os.path.join(local_path, os.path.basename(remote))
    if not os.path.isfile(local):
        log.info(f'downloading KML file to {local_path}')
        with Lock(local_path):
            r = requests.get(remote)
            with open(local, 'wb') as out:
                out.write(r.content)
    try:
        lock = Lock(local_path, soft=True)
        lock.remove()
    except Exception as e:
        raise e
    return local
