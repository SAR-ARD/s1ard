import os
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
        return ['mode', 'aoi_tiles', 'aoi_geometry', 'mindate', 'maxdate', 'acq_mode',
                'work_dir', 'scene_dir', 'rtc_dir', 'tmp_dir', 'wbm_dir', 'measurement',
                'db_file', 'kml_file', 'dem_type', 'gdal_threads', 'log_dir', 'nrb_dir',
                'etad', 'etad_dir', 'product', 'annotation']
    elif section == 'metadata':
        return ['access_url', 'licence', 'doi', 'processing_center']
    else:
        raise RuntimeError(f"unknown section: {section}. Options: 'processing', 'metadata'.")


def get_config(config_file, proc_section='PROCESSING', **kwargs):
    """Returns the content of a `config.ini` file as a dictionary.
    
    Parameters
    ----------
    config_file: str
        Full path to the config file that should be parsed to a dictionary.
    proc_section: str
        Section of the config file that processing parameters should be parsed from. Default is 'PROCESSING'.
    
    Returns
    -------
    out_dict: dict
        Dictionary of the parsed config parameters.
    """
    parser = configparser.ConfigParser(allow_no_value=True,
                                       converters={'_annotation': _parse_annotation,
                                                   '_datetime': _parse_datetime,
                                                   '_tile_list': _parse_tile_list})
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
        proc_sec[k] = v
    
    # set some defaults
    if 'etad' not in proc_sec.keys():
        proc_sec['etad'] = 'False'
        proc_sec['etad_dir'] = 'None'
    for item in ['rtc_dir', 'tmp_dir', 'nrb_dir', 'wbm_dir', 'log_dir']:
        if item not in proc_sec.keys():
            proc_sec[item] = item[:3].upper()
    if 'gdal_threads' not in proc_sec.keys():
        proc_sec['gdal_threads'] = '4'
    if 'dem_type' not in proc_sec.keys():
        proc_sec['dem_type'] = 'Copernicus 30m Global DEM'
    
    # use previous defaults for measurement and annotation if they have not been defined
    if 'measurement' not in proc_sec.keys():
        proc_sec['measurement'] = 'gamma'
    if 'annotation' not in proc_sec.keys():
        if proc_sec['measurement'] == 'gamma':
            proc_sec['annotation'] = 'dm,ei,id,lc,li,np,gs'
        else:
            proc_sec['annotation'] = 'dm,ei,id,lc,li,np,sg'
    
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
            allowed = ['nrb', 'rtc', 'all']
            assert v in allowed, "Parameter '{}': expected to be one of {}; got '{}' instead".format(k, allowed, v)
            v = v.lower()
        if k == 'aoi_tiles':
            if v is not None:
                v = proc_sec.get_tile_list(k)
        if k == 'aoi_geometry':
            if v is not None:
                assert os.path.isfile(v), "Parameter '{}': File {} could not be found".format(k, v)
        if k.endswith('date'):
            v = proc_sec.get_datetime(k)
        if k == 'acq_mode':
            assert v in ['IW', 'EW', 'SM']
        if k == 'work_dir':
            assert os.path.isdir(v), "Parameter '{}': '{}' must be an existing directory".format(k, v)
        dir_ignore = ['work_dir']
        if proc_sec['etad'] == 'False':
            dir_ignore.append('etad_dir')
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
        if k == 'db_file':
            if not any(x in v for x in ['/', '\\']):
                v = os.path.join(proc_sec['work_dir'], v)
        if k == 'gdal_threads':
            v = int(v)
        if k == 'dem_type':
            allowed = ['Copernicus 10m EEA DEM', 'Copernicus 30m Global DEM II',
                       'Copernicus 30m Global DEM', 'GETASSE30']
            assert v in allowed, "Parameter '{}': expected to be one of {}; got '{}' instead".format(k, allowed, v)
        if k == 'etad':
            if v.lower() == 'true':
                v = True
            elif v.lower() == 'false':
                v = False
            else:
                allowed = ['True', 'true', 'False', 'false']
                raise ValueError("Parameter '{}': expected to be one of {}; got '{}' instead".format(k, allowed, v))
        if k == 'product':
            allowed = ['GRD', 'SLC']
            assert v in allowed, "Parameter '{}': expected to be one of {}; got '{}' instead".format(k, allowed, v)
        if k == 'measurement':
            allowed = ['gamma', 'sigma']
            assert v in allowed, "Parameter '{}': expected to be one of {}; got '{}' instead".format(k, allowed, v)
        if k == 'annotation':
            v = proc_sec.get_annotation(k)
        out_dict[k] = v
    
    # METADATA section
    meta_keys = get_keys(section='metadata')
    if 'METADATA' not in parser.keys():
        parser.add_section('METADATA')
    meta_sec = parser['METADATA']
    out_dict['meta'] = {}
    for k, v in meta_sec.items():
        v = _keyval_check(key=k, val=v, allowed_keys=meta_keys)
        # No need to check values. Only requirement is that they're strings, which is configparser's default.
        out_dict['meta'][k] = v
    for key in meta_keys:
        if key not in out_dict['meta'].keys():
            out_dict['meta'][key] = None
    
    return out_dict


def _parse_annotation(s):
    """Custom converter for configparser:
    https://docs.python.org/3/library/configparser.html#customizing-parser-behaviour"""
    if s in ['', 'None']:
        return None
    annotation_list = s.replace(' ', '').split(',')
    allowed = ['dm', 'ei', 'em', 'id', 'lc', 'li', 'np', 'gs', 'sg']
    for annotation in annotation_list:
        if annotation not in allowed:
            msg = "Parameter 'annotation': Error while parsing to list; " \
                  "annotation '{}' is not supported. Allowed keys:\n{}"
            raise ValueError(msg.format(annotation, allowed))
    return annotation_list


def _parse_datetime(s):
    """Custom converter for configparser:
    https://docs.python.org/3/library/configparser.html#customizing-parser-behaviour"""
    return dateutil.parser.parse(s)


def _parse_tile_list(s):
    """Custom converter for configparser:
    https://docs.python.org/3/library/configparser.html#customizing-parser-behaviour"""
    tile_list = s.replace(' ', '').split(',')
    for tile in tile_list:
        if len(tile) != 5:
            raise ValueError("Parameter 'aoi_tiles': Error while parsing MGRS tile IDs to list; tile '{}' is not 5 "
                             "digits long.".format(tile))
        else:
            continue
    return tile_list


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
    Returns a dictionary of additional parameters for :func:`S1_NRB.snap.process` based on processing
    configurations provided by the config file.
    
    Parameters
    ----------
    config: dict
        Dictionary of the parsed config parameters for the current process.
    
    Returns
    -------
    dict
        Dictionary of parameters that can be passed to :func:`S1_NRB.snap.process`
    """
    return {'spacing': {'IW': 10,
                        'SM': 10,
                        'EW': 20}[config['acq_mode']],
            'allow_res_osv': True,
            'dem_resampling_method': 'BILINEAR_INTERPOLATION',
            'img_resampling_method': 'BILINEAR_INTERPOLATION',
            'slc_clean_edges': True,
            'slc_clean_edges_pixels': 4,
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
