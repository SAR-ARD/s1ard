import os
import re
from datetime import timedelta
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


def get_config(config_file, proc_section='PROCESSING', **kwargs):
    """
    Returns the content of a `config.ini` file as a dictionary.
    
    Parameters
    ----------
    config_file: str
        Full path to the config file that should be parsed to a dictionary.
    proc_section: str
        Section of the config file that processing parameters should be parsed from. Default is 'PROCESSING'.
    
    Returns
    -------
    dict
        Dictionary of the parsed config parameters.
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
        parser.add_section(proc_section)
        parser.add_section('METADATA')
    else:
        raise TypeError(f"'config_file' must be of type str or None, was {type(config_file)}")
    out_dict = {}
    
    # PROCESSING section
    allowed_keys = get_keys(section='processing')
    try:
        proc_sec = parser[proc_section]
    except KeyError:
        raise KeyError("Section '{}' does not exist in config file {}".format(proc_section, config_file))
    
    # override config file parameters
    for k, v in kwargs.items():
        proc_sec[k] = v.strip()
    
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
    
    for k, v in proc_sec.items():
        v = _keyval_check(key=k, val=v, allowed_keys=allowed_keys)
        
        if k == 'mode':
            v = proc_sec.get_modes(k)
        if k == 'aoi_tiles':
            if v is not None:
                v = proc_sec.get_tile_list(k)
        if k == 'aoi_geometry':
            if v is not None:
                assert os.path.isfile(v), "Parameter '{}': File {} could not be found".format(k, v)
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
            assert os.path.isdir(v), "Parameter '{}': '{}' must be an existing directory".format(k, v)
        dir_ignore = ['work_dir']
        if proc_sec['etad'] == 'False':
            dir_ignore.append('etad_dir')
        if k == 'scene_dir' and v is None:
            dir_ignore.append(k)
        if k.endswith('_dir') and k not in dir_ignore:
            if any(x in v for x in ['/', '\\']):
                assert os.path.isdir(v), "Parameter '{}': {} is a full path to a non-existing directory".format(k, v)
            else:
                v = os.path.join(proc_sec['work_dir'], v)
        if k.endswith('_file') and not k.startswith('db'):
            if any(x in v for x in ['/', '\\']):
                assert os.path.isfile(v), "Parameter '{}': File {} could not be found".format(k, v)
            else:
                v = os.path.join(proc_sec['work_dir'], v)
                assert os.path.isfile(v), "Parameter '{}': File {} could not be found".format(k, v)
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
            assert v in allowed, "Parameter '{}': expected to be one of {}; got '{}' instead".format(k, allowed, v)
        if k in ['etad', 'date_strict']:
            v = proc_sec.getboolean(k)
        if k == 'product':
            allowed = ['GRD', 'SLC']
            assert v in allowed, "Parameter '{}': expected to be one of {}; got '{}' instead".format(k, allowed, v)
        if k == 'measurement':
            allowed = ['gamma', 'sigma']
            assert v in allowed, "Parameter '{}': expected to be one of {}; got '{}' instead".format(k, allowed, v)
        if k == 'annotation':
            v = proc_sec.get_annotation(k)
        if k == 'snap_gpt_args':
            v = proc_sec.get_list(k)
        if k == 'datatake':
            v = proc_sec.get_list(k)
        out_dict[k] = v
    
    if out_dict['db_file'] is None and out_dict['stac_catalog'] is None:
        raise RuntimeError("Either 'db_file' or 'stac_catalog' has to be defined.")
    if out_dict['db_file'] is not None and out_dict['stac_catalog'] is not None:
        raise RuntimeError("both 'db_file' and 'stac_catalog' have been defined. Please choose only one.")
    if out_dict['stac_catalog'] is not None:
        if out_dict['stac_collections'] is None:
            raise RuntimeError("'stac_collections' must be defined if data is to be searched in a STAC.")
    if out_dict['db_file'] is not None:
        if out_dict['scene_dir'] is None:
            raise RuntimeError("'scene_dir' must be defined if data is to be searched via an SQLite database.")
    
    # METADATA section
    allowed_keys = get_keys(section='metadata')
    if 'METADATA' not in parser.keys():
        parser.add_section('METADATA')
    meta_sec = parser['METADATA']
    out_dict['meta'] = {}
    
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
        out_dict['meta'][k] = v
    for key in allowed_keys:
        if key not in out_dict['meta'].keys():
            out_dict['meta'][key] = None
    
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
                        'EW': 40}[config['acq_mode']],
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
    threads = config['gdal_threads']
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
