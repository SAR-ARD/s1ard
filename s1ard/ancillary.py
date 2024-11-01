import os
import sys
import logging
import requests
import hashlib
import dateutil.parser
from pathlib import Path
from multiformats import multihash
import binascii
from lxml import etree
from textwrap import dedent
from datetime import datetime, timedelta, timezone
from osgeo import gdal
import numpy as np
import spatialist
from spatialist.raster import Raster, rasterize
from spatialist.vector import bbox, intersect, boundary, vectorize
import pyroSAR
from pyroSAR.ancillary import Lock, LockCollection
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


def get_max_ext(geometries, buffer=None, crs=None):
    """
    Gets the maximum extent from a list of geometries.
    
    Parameters
    ----------
    geometries: list[spatialist.vector.Vector]
        List of :class:`~spatialist.vector.Vector` geometries.
    buffer: float or None
        The buffer in units of the geometries' CRS to add to the extent.
    crs: str or int or None
        The target CRS of the extent. If None (default) the extent is
        expressed in the CRS of the input geometries.
    
    Returns
    -------
    max_ext: dict
        The maximum extent of the selected :class:`~spatialist.vector.Vector` geometries including the chosen buffer.
    """
    max_ext = {}
    crs_list = []
    for geo in geometries:
        crs_list.append(f"EPSG:{geo.getProjection('epsg')}")
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
    crs_list = list(set(crs_list))
    if len(crs_list) > 1:
        raise RuntimeError(f'The input geometries are in different CRSs: {crs_list}')
    max_ext = dict(max_ext)
    if buffer is not None:
        max_ext['xmin'] -= buffer
        max_ext['xmax'] += buffer
        max_ext['ymin'] -= buffer
        max_ext['ymax'] += buffer
    if crs is not None:
        with bbox(coordinates=max_ext, crs=crs_list[0]) as geo:
            geo.reproject(projection=crs)
            max_ext = geo.extent
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
    
    logfile = config['processing']['logfile']
    if logfile is not None:
        os.makedirs(os.path.dirname(logfile), exist_ok=True)
        handler = logging.FileHandler(filename=logfile, mode='a')
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
    
    mode                {config['processing']['mode']}
    aoi_tiles           {config['processing']['aoi_tiles']}
    aoi_geometry        {config['processing']['aoi_geometry']}
    scene               {config['processing']['scene']}
    mindate             {config['processing']['mindate'].isoformat()}
    maxdate             {config['processing']['maxdate'].isoformat()}
    date_strict         {config['processing']['date_strict']}
    sensor              {config['processing']['sensor']}
    acq_mode            {config['processing']['acq_mode']}
    product             {config['processing']['product']}
    datatake            {config['processing']['datatake']}
    measurement         {config['processing']['measurement']}
    annotation          {config['processing']['annotation']}
    dem_type            {config['processing']['dem_type']}
    etad                {config['processing']['etad']}
    
    work_dir            {config['processing']['work_dir']}
    sar_dir             {config['processing']['sar_dir']}
    tmp_dir             {config['processing']['tmp_dir']}
    ard_dir             {config['processing']['ard_dir']}
    wbm_dir             {config['processing']['wbm_dir']}
    etad_dir            {config['processing']['etad_dir']}
    scene_dir           {config['processing']['scene_dir']}
    logfile             {config['processing']['logfile']}
    db_file             {config['processing']['db_file']}
    stac_catalog        {config['processing']['stac_catalog']}
    stac_collections    {config['processing']['stac_collections']}
    parquet             {config['processing']['parquet']}
    gdal_threads        {config['processing']['gdal_threads']}
    snap_gpt_args       {config['processing']['snap_gpt_args']}
    
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


def date_to_utc(date, as_datetime=False):
    """
    convert a date object to a UTC date string or datetime object.

    Parameters
    ----------
    date: str or datetime or None
        the date object to convert; timezone-unaware dates are interpreted as UTC.
    as_datetime: bool
        return a datetime object instead of a string?

    Returns
    -------
    str or datetime or None
        the date string or datetime object in UTC time zone
    """
    if date is None:
        return date
    elif isinstance(date, str):
        out = dateutil.parser.parse(date)
    elif isinstance(date, datetime):
        out = date
    else:
        raise TypeError('date must be a string, datetime object or None')
    if out.tzinfo is None:
        out = out.replace(tzinfo=timezone.utc)
    else:
        out = out.astimezone(timezone.utc)
    if not as_datetime:
        out = out.strftime('%Y-%m-%dT%H:%M:%SZ')
    return out


def buffer_time(start, stop, as_datetime=False, **kwargs):
    """
    Time range buffering
    
    Parameters
    ----------
    start: str
        the start time date object to convert; timezone-unaware dates are interpreted as UTC.
    stop: str
        the stop time date object to convert; timezone-unaware dates are interpreted as UTC.
    as_datetime: bool
        return datetime objects instead of strings?
    kwargs
        time arguments passed to :func:`datetime.timedelta`

    Returns
    -------
    tuple[str | datetime]
        the buffered start and stop time as UTC string or datetime object
    """
    td = timedelta(**kwargs)
    start = date_to_utc(start, as_datetime=True) - td
    stop = date_to_utc(stop, as_datetime=True) + td
    if not as_datetime:
        start = start.strftime('%Y-%m-%dT%H:%M:%SZ')
        stop = stop.strftime('%Y-%m-%dT%H:%M:%SZ')
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
    local_path = os.path.join(os.path.expanduser('~'), '.s1ard')
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


def compute_hash(file_path, algorithm='sha256', chunk_size=8192, multihash_encode=True):
    """
    Compute the (multi)hash of a file using the specified algorithm.

    Parameters
    ----------
    file_path: str
        Path to the file.
    algorithm: str
        Hash algorithm to use (default is 'sha256').
    chunk_size: int
        Size of chunks to read from the file in bytes (default is 8192).
    multihash_encode: bool
        Encode the hash according to the
        `multihash specification <https://github.com/multiformats/multihash>`_
        (default is True)?
        The hash generated by `hashlib` will be wrapped using
        :func:`multiformats.multihash.wrap`.

    Returns
    -------
    str
        the hexadecimal hash string of the file.

    See Also
    --------
    :mod:`hashlib`
    :mod:`multiformats.multihash`
    """
    # lookup between hashlib and multihash algorithm names; to be extended if necessary
    algo_lookup = {'sha1': 'sha1',
                   'sha256': 'sha2-256',
                   'sha512': 'sha2-512'}
    if algorithm not in algo_lookup.keys():
        raise ValueError(f'Hash algorithm must be one of {algo_lookup.keys()}')
    hash_func = getattr(hashlib, algorithm)()
    with open(file_path, 'rb') as f:
        while chunk := f.read(chunk_size):
            hash_func.update(chunk)
    if multihash_encode:
        digest = hash_func.digest()
        mh = multihash.wrap(digest, algo_lookup[algorithm])
        return mh.hex()
    else:
        return hash_func.hexdigest()


def datamask(measurement, dm_ras, dm_vec):
    """
    Create data masks for a given image file.
    The created raster data mask does not contain a simple mask of nodata values.
    Rather, a boundary vector geometry containing all valid pixels is created and
    then rasterized. This boundary geometry (single polygon) is saved as `dm_vec`.
    In this case `dm_vec` is returned.
    If the input image only contains nodata values, no raster data mask is created,
    and an empty dummy vector mask is created. In this case the function will return
    `None`.
    
    
    Parameters
    ----------
    measurement: str
        the binary image file
    dm_ras: str
        the name of the raster data mask
    dm_vec: str
        the name of the vector data mask

    Returns
    -------
    str or None
        `dm_vec` if the vector data mask contains a geometry or None otherwise
    """
    
    def mask_from_array(arr, dm_vec, dm_ras, ref):
        """
        
        Parameters
        ----------
        arr: np.ndarray
        dm_vec: str
        dm_ras: str
        ref: spatialist.raster.Raster

        Returns
        -------
        str or None
        """
        # create a dummy vector mask if the mask only contains 0 values
        if len(arr[arr == 1]) == 0:
            Path(dm_vec).touch(exist_ok=False)
            return None
        # vectorize the raster data mask
        with vectorize(target=arr, reference=ref) as vec:
            # compute a valid data boundary geometry (vector data mask)
            with boundary(vec, expression="value=1") as bounds:
                # rasterize the vector data mask
                if not os.path.isfile(dm_ras):
                    rasterize(vectorobject=bounds, reference=ref,
                              outname=dm_ras)
                # write the vector data mask
                bounds.write(outfile=dm_vec)
        return dm_vec
    
    if os.path.isfile(dm_vec) and os.path.isfile(dm_ras):
        return None if os.path.getsize(dm_vec) == 0 else dm_vec
    
    with LockCollection([dm_vec, dm_ras]):
        if not os.path.isfile(dm_vec):
            if not os.path.isfile(dm_ras):
                with Raster(measurement) as ras:
                    arr = ras.array()
                    # create a nodata mask
                    mask = ~np.isnan(arr)
                    del arr
                    out = mask_from_array(arr=mask, dm_vec=dm_vec,
                                          dm_ras=dm_ras, ref=ras)
            else:
                # read the raster data mask
                with Raster(dm_ras) as ras:
                    mask = ras.array()
                    out = mask_from_array(arr=mask, dm_vec=dm_vec,
                                          dm_ras=dm_ras, ref=ras)
                    del mask
        else:
            if os.path.getsize(dm_vec) == 0:
                out = None
            else:
                out = dm_vec
    return out
