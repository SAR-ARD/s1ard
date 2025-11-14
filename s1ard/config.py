import os
import re
import copy
import importlib.resources
from datetime import datetime, timedelta
import configparser
from dateutil.parser import parse as dateparse
from osgeo import gdal
from s1ard.processors.registry import load_processor
from cesard.config import keyval_check, validate_options, validate_value
from typing import Any


def _get_config_metadata(parser, **kwargs):
    # METADATA section
    allowed_keys = get_keys(section='metadata')
    if 'METADATA' not in parser.sections():
        parser.add_section('METADATA')
    meta_sec = parser['METADATA']
    
    # override config file parameters
    for k, v in kwargs.items():
        if k in allowed_keys:
            meta_sec[k] = v.strip()
    
    # set defaults
    if 'format' not in meta_sec.keys():
        meta_sec['format'] = 'OGC, STAC'
    if 'copy_original' not in meta_sec.keys():
        meta_sec['copy_original'] = 'True'
    
    out = {}
    for k, v in meta_sec.items():
        v = keyval_check(key=k, val=v, allowed_keys=allowed_keys)
        if k == 'format':
            v = meta_sec.get_list(k)
        if k == 'copy_original':
            v = meta_sec.getboolean(k)
        out[k] = v
    for key in allowed_keys:
        if key not in out.keys():
            out[key] = None
    return out


def _get_config_processing(parser, **kwargs):
    allowed_keys = get_keys(section='processing')
    try:
        proc_sec = parser['PROCESSING']
    except KeyError:
        msg = "Section 'PROCESSING' does not exist in the config file"
        raise KeyError(msg)
    
    # override config file parameters with additional keyword arguments
    for k, v in kwargs.items():
        if k in allowed_keys:
            proc_sec[k] = v.strip()
    
    # make all relevant paths absolute
    for k in ['work_dir', 'scene_dir', 'scene', 'etad_dir']:
        v = proc_sec[k]
        proc_sec[k] = 'None' if v in ['', 'None'] else os.path.abspath(v)
    
    # set some defaults
    processing_defaults = {
        'sar_dir': 'SAR',
        'tmp_dir': 'TMP',
        'ard_dir': 'ARD',
        'wbm_dir': 'WBM',
        'gdal_threads': '4',
        'dem_type': 'Copernicus 30m Global DEM',
        'date_strict': 'True',
        'datatake': 'None',
        'measurement': 'gamma',
        'annotation': 'dm,ei,id,lc,li,np,ratio',
        'logfile': 'None',
        'parquet': 'None'
    }
    processing_options = {
        'acq_mode': ['IW', 'EW', 'SM'],
        'annotation': ['dm', 'ei', 'em', 'id', 'lc',
                       'ld', 'li', 'np', 'ratio', 'wm'],
        'dem_type': ['Copernicus 10m EEA DEM',
                     'Copernicus 30m Global DEM',
                     'Copernicus 30m Global DEM II',
                     'GETASSE30'],
        'measurement': ['gamma', 'sigma'],
        'mode': ['sar', 'nrb', 'orb'],
        'product': ['GRD', 'SLC'],
        'sensor': ['S1A', 'S1B', 'S1C', 'S1D']}
    
    if 'etad' not in proc_sec.keys():
        proc_sec['etad'] = 'False'
        proc_sec['etad_dir'] = 'None'
    for k, v in processing_defaults.items():
        if k not in proc_sec.keys():
            proc_sec[k] = v
    
    # check completeness of configuration parameters
    missing = []
    exclude = ['aoi_tiles', 'aoi_geometry']
    for key in get_keys(section='processing'):
        if key not in proc_sec.keys() and key not in exclude:
            missing.append(key)
    if len(missing) > 0:
        missing_str = '\n - ' + '\n - '.join(missing)
        raise RuntimeError(f"missing the following parameters:{missing_str}")
    
    out = {}
    for k, v in proc_sec.items():
        # check if key is allowed and convert 'None|none|' strings to None
        v = keyval_check(key=k, val=v, allowed_keys=allowed_keys)
        
        if k in ['annotation', 'aoi_tiles', 'data_take', 'mode', 'stac_collections']:
            v = proc_sec.get_list(k)
        
        validate_value(k, v)
        
        if k == 'mindate' and v is not None:
            v = proc_sec.get_datetime(k)
        if k == 'maxdate' and v is not None:
            date_short = re.search('^[0-9-]{10}$', v) is not None
            v = proc_sec.get_datetime(k)
            if date_short:
                v += timedelta(days=1, microseconds=-1)
        dir_ignore = ['work_dir']
        if proc_sec['etad'] == 'False':
            dir_ignore.append('etad_dir')
        if k == 'scene_dir' and v is None:
            dir_ignore.append(k)
        if k.endswith('_dir') and k not in dir_ignore:
            if os.path.isabs(v):
                msg = f"Parameter '{k}': '{v}' must be an existing directory"
                assert v is not None and os.path.isdir(v), msg
            else:
                v = os.path.join(proc_sec['work_dir'], v)
        if k.endswith('_file') and not k.startswith('db'):
            msg = f"Parameter '{k}': file {v} could not be found"
            if os.path.isabs(v):
                assert os.path.isfile(v), msg
            else:
                v = os.path.join(proc_sec['work_dir'], v)
                assert os.path.isfile(v), msg
        if k in ['db_file', 'logfile'] and v is not None:
            if not os.path.isabs(v):
                v = os.path.join(proc_sec['work_dir'], v)
        if k == 'gdal_threads':
            v = int(v)
        if k in ['etad', 'date_strict']:
            v = proc_sec.getboolean(k)
        
        validate_options(k, v, options=processing_options)
        out[k] = v
    
    # check that a valid scene search option is set
    db_file_set = out['db_file'] is not None
    stac_catalog_set = out['stac_catalog'] is not None
    stac_collections_set = out['stac_collections'] is not None
    parquet_set = out['parquet'] is not None
    
    options_set = sum([db_file_set, stac_catalog_set, parquet_set])
    
    if options_set == 0:
        raise RuntimeError("Please define a scene search option.")
    elif options_set > 1:
        raise RuntimeError("Multiple scene search options have been defined. Please choose only one.")
    
    if stac_catalog_set and not stac_collections_set:
        raise RuntimeError("'stac_collections' must be defined if data is to be searched in a STAC.")
    
    return out


def _parse_datetime(s):
    """Custom converter for configparser:
    https://docs.python.org/3/library/configparser.html#customizing-parser-behaviour"""
    return dateparse(s)


def _parse_list(s):
    """Custom converter for configparser:
    https://docs.python.org/3/library/configparser.html#customizing-parser-behaviour"""
    if s in ['', 'None']:
        return None
    else:
        return [x.strip() for x in s.split(',')]


def gdal_conf(config):
    """
    Stores GDAL configuration options for the current process.

    Parameters
    ----------
    config: dict
        Dictionary of the parsed config parameters for the current process.

    Returns
    -------
    dict
        Dictionary containing GDAL configuration options for the current process.
    """
    threads = config['processing']['gdal_threads']
    threads_before = gdal.GetConfigOption('GDAL_NUM_THREADS')
    if not isinstance(threads, int):
        raise TypeError("'threads' must be of type int")
    if threads == 1:
        multithread = False
    elif threads > 1:
        multithread = True
        gdal.SetConfigOption('GDAL_NUM_THREADS', str(threads))
    else:
        raise ValueError("'threads' must be >= 1")
    
    return {'threads': threads, 'threads_before': threads_before,
            'multithread': multithread}


def get_config(config_file: str | None = None, **kwargs: dict[str, str]) \
        -> dict[str, Any]:
    """
    Returns the content of a `config.ini` file as a dictionary.

    Parameters
    ----------
    config_file:
        Full path to the config file that should be parsed to a dictionary.
    kwargs:
        further keyword arguments overriding configuration found in the config file.

    Returns
    -------
        Dictionary of the parsed config parameters.
        The keys correspond to the config sections in lowercase letters.
    """
    parser = read_config_file(config_file)
    
    kwargs_proc = {k: v for k, v in kwargs.items() if k in get_keys('processing')}
    kwargs_meta = {k: v for k, v in kwargs.items() if k in get_keys('metadata')}
    
    out = {'processing': _get_config_processing(parser, **kwargs_proc),
           'metadata': _get_config_metadata(parser, **kwargs_meta)}
    
    processor_name = out['processing']['processor']
    processor = load_processor(processor_name)
    kwargs_sar = {k: v for k, v in kwargs.items() if k in get_keys(processor_name)}
    out[processor_name] = processor.get_config_section(parser, **kwargs_sar)
    
    return out


def get_keys(section: str) -> list[str]:
    """
    get all allowed configuration keys for a section
    
    Parameters
    ----------
    section:
        the configuration section to get the allowed keys for.
        Either 'processing', 'metadata' or the name of a SAR processor plugin e.g. 'snap'.

    Returns
    -------
        a list of keys
    """
    if section == 'processing':
        return ['acq_mode', 'annotation', 'aoi_geometry', 'aoi_tiles',
                'ard_dir', 'datatake', 'date_strict', 'db_file', 'dem_type',
                'etad', 'etad_dir', 'gdal_threads', 'logfile', 'maxdate',
                'measurement', 'mindate', 'mode', 'parquet', 'processor',
                'product', 'sar_dir', 'scene', 'scene_dir', 'sensor',
                'stac_catalog', 'stac_collections', 'tmp_dir', 'wbm_dir', 'work_dir']
    elif section == 'metadata':
        return ['access_url', 'copy_original', 'doi', 'format', 'licence', 'processing_center']
    else:
        try:
            processor = load_processor(section)
        except ModuleNotFoundError:
            raise RuntimeError(f"unknown section: {section}.")
        try:
            return processor.get_config_keys()
        except AttributeError:
            raise RuntimeError(f"missing function s1ard.{section}.get_config_keys().")


def init(
        target: str,
        source: str | None = None,
        overwrite: bool = False,
        **kwargs
) -> None:
    """
    Initialize a configuration file.

    Parameters
    ----------
    target:
        Path to the target configuration file.
    source:
        Path to the source file to read the configuration from. If not provided,
        a default configuration file within the package will be used.
    overwrite:
        Overwrite an existing file?
    kwargs:
        Additional keyword arguments for overwriting the configuration in `source`.

    Examples
    --------
    Create a file in the current working directory.
    `work_dir` and a scene search option (in this case SQLite via `db_file`)
    must be defined, other configuration is read from the default configuration file.

    >>> from s1ard.config import init
    >>> init(target='config.ini', work_dir='.', db_file='scenes.db')

    """
    if source is None:
        with importlib.resources.path(package='s1ard.resources',
                                      resource='config.ini') as path:
            source = str(path)
    config = get_config(config_file=source, **kwargs)
    write(config=config, target=target, overwrite=overwrite)


def read_config_file(config_file: str | None = None) -> configparser.ConfigParser:
    """
    Reads a configuration file and returns a ConfigParser object
    
    Parameters
    ----------
    config_file: str or None
        the configuration file name. If None, the default configuration file
        within the package will be used.

    Returns
    -------
        the configuration object
    """
    parser = configparser.ConfigParser(allow_no_value=True,
                                       converters={'_datetime': _parse_datetime,
                                                   '_list': _parse_list})
    
    if config_file:
        if not os.path.isfile(config_file):
            raise FileNotFoundError(f"Config file {config_file} does not exist.")
    else:
        with importlib.resources.path(package='s1ard.resources',
                                      resource='config.ini') as path:
            config_file = str(path)
    
    parser.read(config_file)
    return parser


def write(config, target, overwrite=False, **kwargs):
    """
    Write configuration options to a config file.
    
    Parameters
    ----------
    config: dict
        the configuration as returned by :func:`get_config`
    target: str
        the name of the output file
    overwrite: bool
        overwrite an existing file if it exists?
    kwargs
        further keyword arguments overriding configuration found in `config`.

    Returns
    -------

    """
    if os.path.isfile(target) and not overwrite:
        raise RuntimeError("target already exists")
    
    def to_string(item):
        """
        
        Parameters
        ----------
        item: dict or List or str

        Returns
        -------
        str or dict
        """
        if isinstance(item, dict):
            return {k: to_string(v) for k, v in item.items()}
        elif isinstance(item, list):
            return ', '.join([to_string(x) for x in item])
        elif isinstance(item, datetime):
            return item.strftime('%Y-%m-%d %H:%M:%S')
        else:
            return str(item)
    
    processor_name = config['processing']['processor']
    processor = load_processor(processor_name)
    
    config = copy.deepcopy(config)
    keys_processing = get_keys('processing')
    keys_meta = get_keys('metadata')
    keys_proc = processor.get_config_keys()
    for k, v in kwargs.items():
        if k in keys_processing:
            config['processing'][k] = v
        elif k in keys_meta:
            config['metadata'][k] = v
        elif k in keys_proc:
            config[processor_name][k] = v
        else:
            raise KeyError("Parameter '{}' is not supported".format(k))
    keys_path_relative = ['sar_dir', 'tmp_dir', 'ard_dir', 'wbm_dir', 'db_file']
    work_dir = config['processing']['work_dir']
    for k in keys_path_relative:
        v = config['processing'][k]
        if v is not None and work_dir in v:
            config['processing'][k] = v.replace(work_dir, '').strip('/\\')
    config['metadata'] = to_string(config['metadata'])
    config['processing'] = to_string(config['processing'])
    config_proc_str = processor.config_to_string(config[processor_name])
    config[processor_name] = config_proc_str
    parser = configparser.ConfigParser()
    parser['METADATA'] = config['metadata']
    parser['PROCESSING'] = config['processing']
    parser[processor_name.upper()] = config[processor_name]
    with open(target, 'w') as configfile:
        parser.write(configfile)
