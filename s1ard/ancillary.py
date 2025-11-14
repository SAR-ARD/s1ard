import os
import sys
import logging
from datetime import datetime
from osgeo import gdal
import spatialist
import pyroSAR
import s1ard
from s1ard.processors.registry import load_processor

log = logging.getLogger('s1ard')


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
    
    log_format = "[%(asctime)s] [%(levelname)7s] %(message)s"
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


def _log_process_config(logger, config):
    """
    Adds a header to the logfile, which includes information about
    the current processing configuration.
    
    Parameters
    ----------
    logger: logging.Logger
        The logger to which the header is added to.
    config: dict
        Dictionary of the parsed config parameters for the current process.
    """
    sw_versions = {
        's1ard': s1ard.__version__,
        'python': sys.version,
        'python-pyroSAR': pyroSAR.__version__,
        'python-spatialist': spatialist.__version__,
        'python-GDAL': gdal.__version__}
    
    processor_name = config['processing']['processor']
    processor = load_processor(processor_name)
    sw_versions.update(processor.version_dict())
    
    max_len_sw = len(max(sw_versions.keys(), key=len))
    max_len_main = len(max(config['processing'].keys(), key=len))
    max_len_meta = len(max(config['metadata'].keys(), key=len))
    max_len_proc = len(max(config[processor_name].keys(), key=len))
    max_len = max(max_len_sw, max_len_main, max_len_meta, max_len_proc) + 4
    
    lines = []
    lines.append('=' * 100)
    for section in ['PROCESSING', processor_name.upper(), 'METADATA']:
        lines.append(f'{section}')
        for k, v in config[section.lower()].items():
            if k == 'dem_prepare_mode':
                continue
            if isinstance(v, datetime):
                val = v.isoformat()
            else:
                val = v
            lines.append(f"{k: <{max_len}}{val}")
        lines.append('=' * 100)
    lines.append('SOFTWARE')
    for k, v in sw_versions.items():
        lines.append(f"{k: <{max_len}}{v}")
    lines.append('=' * 100)
    header = '\n'.join(lines)
    logger.info(header)
