import os
import configparser
from datetime import datetime
from osgeo import gdal


def get_config(config_file, proc_section='PROCESSING'):
    """Returns the content of a `config.ini` file as a dictionary.
    
    Parameters
    ----------
    config_file: str
        Full path to the config file that should be parsed to a dictionary.
    proc_section: str, optional
        Section of the config file that processing parameters should be parsed from. Default is 'PROCESSING'.
    
    Returns
    -------
    out_dict: dict
        Dictionary of the parsed config parameters.
    """
    if not os.path.isfile(config_file):
        raise FileNotFoundError("Config file {} does not exist.".format(config_file))
    
    parser = configparser.ConfigParser(allow_no_value=True, converters={'_datetime': _parse_datetime,
                                                                        '_tile_list': _parse_tile_list})
    parser.read(config_file)
    out_dict = {}
    
    # PROCESSING section
    allowed_keys = ['mode', 'aoi_tiles', 'aoi_geometry', 'mindate', 'maxdate', 'acq_mode',
                    'work_dir', 'scene_dir', 'rtc_dir', 'tmp_dir', 'dem_dir', 'wbm_dir',
                    'db_file', 'kml_file', 'dem_type', 'gdal_threads', 'log_dir', 'nrb_dir',
                    'etad', 'etad_dir', 'product']
    try:
        proc_sec = parser[proc_section]
    except KeyError:
        raise KeyError("Section '{}' does not exist in config file {}".format(proc_section, config_file))
    
    for k, v in proc_sec.items():
        v = _keyval_check(key=k, val=v, allowed_keys=allowed_keys)
        
        if k == 'mode':
            allowed = ['nrb', 'snap', 'all']
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
                os.makedirs(v, exist_ok=True)
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
        out_dict[k] = v
    
    assert any([out_dict[k] is not None for k in ['aoi_tiles', 'aoi_geometry']])
    
    # METADATA section
    meta_keys = ['access_url', 'licence', 'doi', 'processing_center']
    try:
        meta_sec = parser['METADATA']
        out_dict['meta'] = {}
        for k, v in meta_sec.items():
            v = _keyval_check(key=k, val=v, allowed_keys=meta_keys)
            # No need to check values. Only requirement is that they're strings, which is configparser's default.
            out_dict['meta'][k] = v
    except KeyError:
        # Use None for all relevant fields if the metadata section doesn't exist.
        out_dict['meta'] = dict([(k, None) for k in meta_keys])
    
    return out_dict


def _parse_datetime(s):
    """Custom converter for configparser:
    https://docs.python.org/3/library/configparser.html#customizing-parser-behaviour"""
    if 'T' in s:
        try:
            return datetime.strptime(s, '%Y-%m-%dT%H:%M:%S')
        except ValueError as e:
            raise Exception("Parameters 'mindate/maxdate': Could not parse '{}' with datetime format "
                            "'%Y-%m-%dT%H:%M:%S'".format(s)) from e
    else:
        try:
            return datetime.strptime(s, '%Y-%m-%d')
        except ValueError as e:
            raise Exception("Parameters 'mindate/maxdate': Could not parse '{}' with datetime format "
                            "'%Y-%m-%d'".format(s)) from e


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


def geocode_conf(config):
    """
    Returns a dictionary of additional parameters for :func:`pyroSAR.snap.util.geocode` based on processing
    configurations provided by the config file.
    
    Parameters
    ----------
    config: dict
        Dictionary of the parsed config parameters for the current process.
    
    Returns
    -------
    dict
        Dictionary of parameters that can be passed to :func:`pyroSAR.snap.util.geocode`
    """
    return {'spacing': {'IW': 10,
                        'SM': 10,
                        'EW': 20}[config['acq_mode']],
            'scaling': 'linear',
            'groupsize': 1,
            'allow_RES_OSV': True,
            'alignToStandardGrid': True,
            'export_extra': ['localIncidenceAngle', 'incidenceAngleFromEllipsoid',
                             'scatteringArea', 'layoverShadowMask', 'gammaSigmaRatio'],
            'refarea': ['sigma0', 'gamma0'],
            'externalDEMApplyEGM': False,
            'demResamplingMethod': 'BILINEAR_INTERPOLATION',
            'clean_edges': True,
            'clean_edges_npixels': 4,
            'test': False,
            'cleanup': True}


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
