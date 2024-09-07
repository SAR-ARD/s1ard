import os
import re
import copy
import importlib.resources
from datetime import datetime, timedelta
import configparser
import dateutil.parser
from osgeo import gdal


def get_keys(section):
    """
    get all allowed configuration keys
    
    Parameters
    ----------
    section: {'processing', 'metadata'}
        the configuration section to get the allowed keys for.

    Returns
    -------
    list[str]
        a list of keys
    """
    if section == 'processing':
        return ['mode', 'aoi_tiles', 'aoi_geometry', 'mindate', 'maxdate', 'acq_mode', 'datatake',
                'work_dir', 'scene_dir', 'sar_dir', 'tmp_dir', 'wbm_dir', 'measurement',
                'db_file', 'dem_type', 'gdal_threads', 'logfile', 'ard_dir',
                'etad', 'etad_dir', 'product', 'annotation', 'stac_catalog', 'stac_collections',
                'sensor', 'date_strict', 'snap_gpt_args', 'scene']
    elif section == 'metadata':
        return ['format', 'copy_original', 'access_url', 'licence', 'doi', 'processing_center']
    else:
        raise RuntimeError(f"unknown section: {section}. Options: 'processing', 'metadata'.")


def get_config(config_file=None, **kwargs):
    """
    Returns the content of a `config.ini` file as a dictionary.
    
    Parameters
    ----------
    config_file: str or None
        Full path to the config file that should be parsed to a dictionary.
    kwargs
        further keyword arguments overriding configuration found in the config file.
    
    Returns
    -------
    dict
        Dictionary of the parsed config parameters.
        The keys correspond to the config sections in lowercase letters.
    """
    parser = configparser.ConfigParser(allow_no_value=True,
                                       converters={'_annotation': _parse_annotation,
                                                   '_datetime': _parse_datetime,
                                                   '_modes': _parse_modes,
                                                   '_stac_collections': _parse_list,
                                                   '_tile_list': _parse_tile_list,
                                                   '_list': _parse_list})
    if isinstance(config_file, str):
        if not os.path.isfile(config_file):
            raise FileNotFoundError("Config file {} does not exist.".format(config_file))
        parser.read(config_file)
    elif config_file is None:
        with importlib.resources.path('s1ard.resources', 'config.ini') as path:
            config_file = str(path)
        parser.read(config_file)
    else:
        raise TypeError(f"'config_file' must be of type str or None, was {type(config_file)}")
    out_dict = {'processing': {},
                'metadata': {}}
    
    # PROCESSING section
    allowed_keys = get_keys(section='processing')
    try:
        proc_sec = parser['PROCESSING']
    except KeyError:
        msg = "Section '{}' does not exist in config file {}"
        raise KeyError(msg.format('PROCESSING', config_file))
    
    # override config file parameters with additional keyword arguments
    for k, v in kwargs.items():
        if k in allowed_keys:
            proc_sec[k] = v.strip()
    
    # make all relevant paths absolute
    for k in ['work_dir', 'scene_dir', 'scene', 'logfile', 'etad_dir']:
        v = proc_sec[k]
        proc_sec[k] = 'None' if v in ['', 'None'] else os.path.abspath(v)
    
    # set some defaults
    if 'etad' not in proc_sec.keys():
        proc_sec['etad'] = 'False'
        proc_sec['etad_dir'] = 'None'
    for item in ['sar_dir', 'tmp_dir', 'ard_dir', 'wbm_dir']:
        if item not in proc_sec.keys():
            proc_sec[item] = item[:3].upper()
    if 'gdal_threads' not in proc_sec.keys():
        proc_sec['gdal_threads'] = '4'
    if 'dem_type' not in proc_sec.keys():
        proc_sec['dem_type'] = 'Copernicus 30m Global DEM'
    if 'date_strict' not in proc_sec.keys():
        proc_sec['date_strict'] = 'True'
    if 'snap_gpt_args' not in proc_sec.keys():
        proc_sec['snap_gpt_args'] = 'None'
    if 'datatake' not in proc_sec.keys():
        proc_sec['datatake'] = 'None'
    # use previous defaults for measurement and annotation if they have not been defined
    if 'measurement' not in proc_sec.keys():
        proc_sec['measurement'] = 'gamma'
    if 'annotation' not in proc_sec.keys():
        proc_sec['annotation'] = 'dm,ei,id,lc,li,np,ratio'
    if 'logfile' not in proc_sec.keys():
        proc_sec['logfile'] = 'None'
    
    # check completeness of configuration parameters
    missing = []
    exclude = ['aoi_tiles', 'aoi_geometry']
    for key in get_keys(section='processing'):
        if key not in proc_sec.keys() and key not in exclude:
            missing.append(key)
    if len(missing) > 0:
        missing_str = '\n - ' + '\n - '.join(missing)
        raise RuntimeError(f"missing the following parameters:{missing_str}")
    
    # convert values to Python objects and validate them
    for k, v in proc_sec.items():
        # check if key is allowed and convert 'None|none|' strings to None
        v = _keyval_check(key=k, val=v, allowed_keys=allowed_keys)
        
        if k == 'mode':
            v = proc_sec.get_modes(k)
        if k == 'aoi_tiles':
            if v is not None:
                v = proc_sec.get_tile_list(k)
        if k == 'aoi_geometry':
            if v is not None:
                msg = f"Parameter '{k}': File {v} could not be found"
                assert os.path.isfile(v), msg
        if k == 'mindate':
            v = proc_sec.get_datetime(k)
        if k == 'maxdate':
            date_short = re.search('^[0-9-]{10}$', v) is not None
            v = proc_sec.get_datetime(k)
            if date_short:
                v += timedelta(days=1, microseconds=-1)
        if k == 'sensor':
            assert v in ['S1A', 'S1B']
        if k == 'acq_mode':
            assert v in ['IW', 'EW', 'SM']
        if k == 'work_dir':
            msg = f"Parameter '{k}': '{v}' must be an existing writable directory"
            assert v is not None and os.path.isdir(v) and os.access(v, os.W_OK), msg
        dir_ignore = ['work_dir']
        if proc_sec['etad'] == 'False':
            dir_ignore.append('etad_dir')
        if k == 'scene_dir' and v is None:
            dir_ignore.append(k)
        if k.endswith('_dir') and k not in dir_ignore:
            if any(x in v for x in ['/', '\\']):
                msg = f"Parameter '{k}': '{v}' must be an existing directory"
                assert v is not None and os.path.isdir(v), msg
            else:
                v = os.path.join(proc_sec['work_dir'], v)
        if k.endswith('_file') and not k.startswith('db'):
            if any(x in v for x in ['/', '\\']):
                msg = f"Parameter '{k}': file {v} could not be found"
                assert os.path.isfile(v), msg
            else:
                v = os.path.join(proc_sec['work_dir'], v)
                msg = f"Parameter '{k}': file {v} could not be found"
                assert os.path.isfile(v), msg
        if k in ['db_file', 'logfile'] and v is not None:
            if not any(x in v for x in ['/', '\\']):
                v = os.path.join(proc_sec['work_dir'], v)
        if k == 'stac_collections':
            v = proc_sec.get_stac_collections(k)
        if k == 'gdal_threads':
            v = int(v)
        if k == 'dem_type':
            allowed = ['Copernicus 10m EEA DEM', 'Copernicus 30m Global DEM II',
                       'Copernicus 30m Global DEM', 'GETASSE30']
            msg = "Parameter '{}': expected to be one of {}; got '{}' instead"
            assert v in allowed, msg.format(k, allowed, v)
        if k in ['etad', 'date_strict']:
            v = proc_sec.getboolean(k)
        if k == 'product':
            allowed = ['GRD', 'SLC']
            msg = "Parameter '{}': expected to be one of {}; got '{}' instead"
            assert v in allowed, msg.format(k, allowed, v)
        if k == 'measurement':
            allowed = ['gamma', 'sigma']
            msg = "Parameter '{}': expected to be one of {}; got '{}' instead"
            assert v in allowed, msg.format(k, allowed, v)
        if k == 'annotation':
            v = proc_sec.get_annotation(k)
        if k == 'snap_gpt_args':
            v = proc_sec.get_list(k)
        if k == 'datatake':
            v = proc_sec.get_list(k)
        out_dict['processing'][k] = v
    
    # check that a valid scene search option is set
    db_file_set = out_dict['processing']['db_file'] is not None
    stac_catalog_set = out_dict['processing']['stac_catalog'] is not None
    stac_collections_set = out_dict['processing']['stac_collections'] is not None
    
    if not db_file_set and not stac_catalog_set:
        raise RuntimeError("Either 'db_file' or 'stac_catalog' has to be defined.")
    if db_file_set and stac_catalog_set:
        raise RuntimeError("both 'db_file' and 'stac_catalog' have been defined. Please choose only one.")
    if stac_catalog_set and not stac_collections_set:
        raise RuntimeError("'stac_collections' must be defined if data is to be searched in a STAC.")
    
    # METADATA section
    allowed_keys = get_keys(section='metadata')
    if 'METADATA' not in parser.keys():
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
    
    for k, v in meta_sec.items():
        v = _keyval_check(key=k, val=v, allowed_keys=allowed_keys)
        if k == 'format':
            v = meta_sec.get_list(k)
        if k == 'copy_original':
            v = meta_sec.getboolean(k)
        out_dict['metadata'][k] = v
    for key in allowed_keys:
        if key not in out_dict['metadata'].keys():
            out_dict['metadata'][key] = None
    
    return out_dict


def _parse_annotation(s):
    """Custom converter for configparser:
    https://docs.python.org/3/library/configparser.html#customizing-parser-behaviour"""
    annotation_list = _parse_list(s)
    if annotation_list is not None:
        allowed = ['dm', 'ei', 'em', 'id', 'lc', 'ld', 'li', 'np', 'ratio', 'wm']
        for layer in annotation_list:
            if layer not in allowed:
                msg = "Parameter 'annotation': Error while parsing to list; " \
                      "layer '{}' is not supported. Allowed keys:\n{}"
                raise ValueError(msg.format(layer, allowed))
    return annotation_list


def _parse_datetime(s):
    """Custom converter for configparser:
    https://docs.python.org/3/library/configparser.html#customizing-parser-behaviour"""
    return dateutil.parser.parse(s)


def _parse_modes(s):
    """Custom converter for configparser:
    https://docs.python.org/3/library/configparser.html#customizing-parser-behaviour"""
    mode_list = _parse_list(s)
    allowed = ['sar', 'nrb', 'orb']
    for mode in mode_list:
        if mode not in allowed:
            msg = "Parameter 'annotation': Error while parsing to list; " \
                  "mode '{}' is not supported. Allowed keys:\n{}"
            raise ValueError(msg.format(mode, allowed))
    return mode_list


def _parse_tile_list(s):
    """Custom converter for configparser:
    https://docs.python.org/3/library/configparser.html#customizing-parser-behaviour"""
    tile_list = _parse_list(s)
    if tile_list is not None:
        for tile in tile_list:
            if len(tile) != 5:
                raise ValueError("Parameter 'aoi_tiles': Error while parsing "
                                 "MGRS tile IDs to list; tile '{}' is not 5 "
                                 "digits long.".format(tile))
            else:
                continue
    return tile_list


def _parse_list(s):
    """Custom converter for configparser:
    https://docs.python.org/3/library/configparser.html#customizing-parser-behaviour"""
    if s in ['', 'None']:
        return None
    else:
        return [x.strip() for x in s.split(',')]


def _keyval_check(key, val, allowed_keys):
    """Helper function to check and clean up key,value pairs while parsing a config file."""
    if key not in allowed_keys:
        raise ValueError("Parameter '{}' is not allowed; should be one of {}".format(key, allowed_keys))
    
    val = val.replace('"', '').replace("'", "")
    if val in ['None', 'none', '']:
        val = None
    
    return val


def snap_conf(config):
    """
    Returns a dictionary of additional parameters for :func:`s1ard.snap.process` based on processing
    configurations provided by the config file.
    
    Parameters
    ----------
    config: dict
        Dictionary of the parsed config parameters for the current process.
    
    Returns
    -------
    dict
        Dictionary of parameters that can be passed to :func:`s1ard.snap.process`
    """
    return {'spacing': {'IW': 10,
                        'SM': 10,
                        'EW': 40}[config['processing']['acq_mode']],
            'allow_res_osv': True,
            'dem_resampling_method': 'BILINEAR_INTERPOLATION',
            'img_resampling_method': 'BILINEAR_INTERPOLATION',
            'clean_edges': True,
            'clean_edges_pixels': 4,
            'cleanup': True
            }


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
        overwrite existing file if it exists?
    kwargs
        further keyword arguments overriding configuration found in the config file.

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
    
    config = copy.deepcopy(config)
    keys_processing = get_keys('processing')
    keys_meta = get_keys('metadata')
    for k, v in kwargs.items():
        if k in keys_processing:
            config['processing'][k] = v
        elif k in keys_meta:
            config['metadata'][k] = v
        else:
            raise KeyError("Parameter '{}' is not supported".format(k))
    keys_path_relative = ['sar_dir', 'tmp_dir', 'ard_dir', 'wbm_dir', 'db_file']
    work_dir = config['processing']['work_dir']
    for k in keys_path_relative:
        v = config['processing'][k]
        if v is not None and work_dir in v:
            config['processing'][k] = v.replace(work_dir, '').strip('/\\')
    config = to_string(config)
    parser = configparser.ConfigParser()
    parser['METADATA'] = config['metadata']
    parser['PROCESSING'] = config['processing']
    with open(target, 'w') as configfile:
        parser.write(configfile)
