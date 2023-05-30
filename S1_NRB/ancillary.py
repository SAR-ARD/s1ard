import os
import sys
import logging
import binascii
from lxml import etree
from datetime import datetime, timedelta
from osgeo import gdal
import spatialist
import pyroSAR
from pyroSAR import examine, identify_many
import S1_NRB


def check_acquisition_completeness(scenes, archive):
    """
    Check presence of neighboring acquisitions.
    Check that for each scene a predecessor and successor can be queried
    from the database unless the scene is at the start or end of the data take.
    This ensures that no scene that could be covering an area of interest is missed
    during processing.
    
    Parameters
    ----------
    scenes: list[pyroSAR.drivers.ID]
        a list of scenes
    archive: pyroSAR.drivers.Archive
        an open scene archive connection

    Returns
    -------
    
    Raises
    ------
    RuntimeError
    """
    messages = []
    for scene in scenes:
        print(scene.scene)
        slice = scene.meta['sliceNumber']
        n_slices = scene.meta['totalSlices']
        groupsize = 3
        has_successor = True
        has_predecessor = True
        if slice == 1:  # first slice in the data take
            groupsize -= 1
            has_predecessor = False
        if slice == n_slices:  # last slice in the data take
            groupsize -= 1
            has_successor = False
        f = '%Y%m%dT%H%M%S'
        td = timedelta(seconds=2)
        start = datetime.strptime(scene.start, f) - td
        start = datetime.strftime(start, f)
        stop = datetime.strptime(scene.stop, f) + td
        stop = datetime.strftime(stop, f)
        # Do another database selection to get the scene in question as well as its potential
        # predecessor and successor by adding an acquisition time buffer of two seconds.
        group = archive.select(sensor=scene.sensor,
                               product=scene.product,
                               acquisition_mode=scene.acquisition_mode,
                               mindate=start,
                               maxdate=stop,
                               date_strict=False)
        group = identify_many(group)
        # if the number of selected scenes is lower than the expected group size,
        # check whether the predecessor, the successor or both are missing.
        if len(group) < groupsize:
            start_min = min([x.start for x in group])
            stop_max = max([x.stop for x in group])
            missing = []
            if start_min > start and has_predecessor:
                missing.append('predecessor')
            if stop_max < stop and has_successor:
                missing.append('successor')
            base = os.path.basename(scene.scene)
            messages.append(f'{" and ".join(missing)} acquisition for scene {base}')
    if len(messages) != 0:
        text = '\n - '.join(messages)
        raise RuntimeError(f'missing the following scenes:\n - {text}')


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
    if 9780 % spacing != 0:
        raise RuntimeError(f'target spacing of {spacing} m does not align with MGRS tile overlap of 9780 m.')


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
        Set pyroSAR logging level to DEBUG?
    
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
    
    mode                {config['mode']}
    aoi_tiles           {config['aoi_tiles']}
    aoi_geometry        {config['aoi_geometry']}
    mindate             {config['mindate'].isoformat()}
    maxdate             {config['maxdate'].isoformat()}
    acq_mode            {config['acq_mode']}
    product             {config['product']}
    measurement         {config.get('measurement')}
    annotation          {config.get('annotation')}
    dem_type            {config.get('dem_type')}
    etad                {config.get('etad')}
    
    work_dir            {config['work_dir']}
    scene_dir           {config['scene_dir']}
    rtc_dir             {config['rtc_dir']}
    tmp_dir             {config['tmp_dir']}
    nrb_dir             {config['nrb_dir']}
    wbm_dir             {config['wbm_dir']}
    log_dir             {config['log_dir']}
    etad_dir            {config['etad_dir']}
    db_file             {config['db_file']}
    kml_file            {config['kml_file']}
    gdal_threads        {config.get('gdal_threads')}
    snap_gpt_args       {config['snap_gpt_args']}
    
    ====================================================================================================================
    SOFTWARE
    
    S1_NRB              {S1_NRB.__version__}
    snap-core           {snap_core}
    snap-s1tbx          {snap_s1tbx}
    python              {sys.version}
    python-pyroSAR      {pyroSAR.__version__}
    python-spatialist   {spatialist.__version__}
    python-GDAL         {gdal.__version__}
    
    ====================================================================================================================
    """
    print(header)
    logger.info(header)


def log(handler, mode, proc_step, scenes, msg):
    """
    Format and handle log messages during processing.
    
    Parameters
    ----------
    handler: logging.Logger
        The log handler for the current process.
    mode: {'info', 'warning', 'exception'}
        Calls the respective logging helper function. E.g., ``handler.info()``.
    proc_step: str
        The processing step for which the message is logged.
    scenes: str or list[str]
        Scenes that are currently being processed.
    msg: Any
        The message that should be logged.
    """
    proc_step = proc_step.zfill(7).replace('0', ' ')
    message = '[{proc_step}] -- {scenes} -- {msg}'
    message = message.format(proc_step=proc_step, scenes=scenes, msg=msg)
    if mode == 'info':
        handler.info(message)
    elif mode == 'warning':
        handler.warning(message)
    elif mode == 'exception':
        handler.exception(message)
    else:
        raise RuntimeError('log mode {} is not supported'.format(mode))


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
