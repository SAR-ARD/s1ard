import os
import re
import shutil
import zipfile
import math
from typing import Union
from statistics import mean
from lxml import etree
from dateutil.parser import parse as dateparse
from datetime import timezone
import numpy as np
from spatialist import Raster
from spatialist.auxil import gdalwarp
from spatialist.ancillary import finder, dissolve
from spatialist.raster import rasterize
from pyroSAR.drivers import ID
from osgeo import gdal
import s1ard
from s1ard.metadata.mapping import (RES_MAP_SLC, RES_MAP_GRD,
                                    ENL_MAP_GRD, OSV_MAP, SLC_ACC_MAP,
                                    URL as URL_PACKAGE)
from s1ard.processors.registry import load_processor
from cesard.ancillary import get_tmp_name
from cesard.metadata.extract import (calc_enl, calc_performance_estimates,
                                     geometry_from_vec, vec_from_srccoords)
from cesard.metadata.mapping import DEM_MAP, LERC_ERR_THRES, URL as URL_BASE

gdal.UseExceptions()

URL = URL_BASE
URL.update(URL_PACKAGE)


def append_wind_norm(meta, wm_ref_speed, wm_ref_direction):
    """
    Update a metadata dictionary with wind model information

    Parameters
    ----------
    meta: dict
        metadata extracted by :func:`meta_dict`
    wm_ref_speed: List[str]
        List of paths pointing to the wind model reference speed files.
    wm_ref_direction: List[str]
        List of paths pointing to the wind model reference direction files.

    Returns
    -------

    """
    if wm_ref_speed is not None and wm_ref_direction is not None:
        wm_ref_mean_speed, wm_ref_mean_dir = calc_wm_ref_stats(wm_ref_speed=wm_ref_speed,
                                                               wm_ref_direction=wm_ref_direction,
                                                               epsg=meta['prod']['crsEPSG'],
                                                               bounds=meta['prod']['geom_stac_bbox_native'])
        meta['prod']['windNormBackscatterMeasurement'] = 'sigma0'
        meta['prod']['windNormBackscatterConvention'] = 'intensity ratio'
        meta['prod']['windNormReferenceDirection'] = wm_ref_mean_dir
        meta['prod']['windNormReferenceModel'] = URL['windNormReferenceModel']
        meta['prod']['windNormReferenceSpeed'] = wm_ref_mean_speed
        meta['prod']['windNormReferenceType'] = 'sigma0-ref'
    else:
        meta['prod']['windNormBackscatterMeasurement'] = None
        meta['prod']['windNormBackscatterConvention'] = None
        meta['prod']['windNormReferenceDirection'] = None
        meta['prod']['windNormReferenceModel'] = None
        meta['prod']['windNormReferenceSpeed'] = None
        meta['prod']['windNormReferenceType'] = None


def calc_geolocation_accuracy(swath_identifier, ei_tif, etad, decimals=2):
    """
    Calculates the radial root mean square error, which is a target requirement of the CARD4L NRB specification
    (Item 4.3). For more information see: https://s1ard.readthedocs.io/en/latest/general/geoaccuracy.html.
    Currently only the Copernicus DEM is supported.

    Parameters
    ----------
    swath_identifier: str
        Swath identifier dependent on acquisition mode.
    ei_tif: str
        Path to the annotation GeoTIFF layer 'Ellipsoidal Incident Angle' of the current product.
    etad: bool
        Was the ETAD correction applied?
    decimals: int, optional
        Number of decimal places to round the calculated rRMSE value to. Default is 2.

    Returns
    -------
    rmse_planar: float or None
        The calculated rRMSE value rounded to two decimal places or None if a DEM other than Copernicus is used.
    """
    if etad:
        # https://sentinel.esa.int/nl/web/sentinel/missions/sentinel-1/data-products/etad-dataset
        slc_acc = {'ALE': {'rg': 0,
                           'az': 0},
                   '1sigma': {'rg': 0.2,
                              'az': 0.1}}
    else:
        swath_id = 'SM' if re.search('S[1-6]', swath_identifier) is not None else swath_identifier
        slc_acc = SLC_ACC_MAP[swath_id]
    
    # min/max ellipsoidal incidence angle
    with Raster(ei_tif) as ras:
        stats = ras.allstats(approximate=False)
        ei_min = stats[0]['min']
    
    if ei_min == 0:
        raise RuntimeError(f'minimum ellipsoid incidence angle cannot be 0\n'
                           f'(file: {ei_tif})')
    
    # Remove generated '.aux.xml' file
    aux = finder(os.path.dirname(ei_tif), ['.tif.aux.xml$'], regex=True, recursive=False)
    for file in aux:
        os.remove(file)
    
    # COP-DEM global mean accuracy (LE68) based on LE90 under assumption of gaussian distribution:
    copdem_glob_1sigma_le68 = 1.56
    rmse_dem_planar = copdem_glob_1sigma_le68 / math.tan(math.radians(ei_min))
    
    # Calculation of RMSE_planar
    rmse_rg = math.sqrt(slc_acc['ALE']['rg'] ** 2 + slc_acc['1sigma']['rg'] ** 2)
    rmse_az = math.sqrt(slc_acc['ALE']['az'] ** 2 + slc_acc['1sigma']['az'] ** 2)
    
    rmse_planar = math.sqrt((rmse_rg / math.sin(math.radians(ei_min))) ** 2 +
                            rmse_az ** 2 +
                            rmse_dem_planar ** 2)
    
    return round(rmse_planar, decimals)


def calc_pslr_islr(annotation_dict, decimals=2):
    """
    Extracts all values for Peak Side Lobe Ratio (PSLR) and Integrated Side Lobe Ratio (ISLR) from the annotation
    metadata of a scene and calculates the mean value for all swaths.

    Parameters
    ----------
    annotation_dict: dict
        A dictionary of annotation files in the form: {'swath ID':`lxml.etree._Element` object}
    decimals: int, optional
        Number of decimal places to round the calculated values to. Default is 2.

    Returns
    -------
    tuple[float]
        a tuple with the following values:

        - pslr: Mean PSLR value for all swaths of the scene.
        - islr: Mean ISLR value for all swaths of the scene.
    """
    swaths = list(annotation_dict.keys())
    pslr_dict = find_in_annotation(annotation_dict=annotation_dict, pattern='.//crossCorrelationPslr', out_type='float')
    islr_dict = find_in_annotation(annotation_dict=annotation_dict, pattern='.//crossCorrelationIslr', out_type='float')
    
    # Mean values per swath
    pslr_mean = {}
    islr_mean = {}
    for swath in swaths:
        pslr_mean[swath] = np.nanmean(pslr_dict[swath])
        islr_mean[swath] = np.nanmean(islr_dict[swath])
    
    # Mean value for all swaths
    pslr = np.round(np.nanmean(list(pslr_mean.values())), decimals)
    islr = np.round(np.nanmean(list(islr_mean.values())), decimals)
    return pslr, islr


def calc_wm_ref_stats(wm_ref_speed, wm_ref_direction, epsg, bounds, resolution=915):
    """
    Calculates the mean wind model reference speed and direction for the wind model annotation layer.

    Parameters
    ----------
    wm_ref_speed: list[str]
        List of paths pointing to the wind model reference speed files.
    wm_ref_speed: list[str]
        List of paths pointing to the wind model reference direction files.
    epsg: int
        The EPSG code of the current MGRS tile.
    bounds: list[float]
        The bounds of the current MGRS tile.
    resolution: int, optional
        The resolution of the wind model reference files in meters. Default is 915.

    Returns
    -------
    tuple[float]
        a tuple with the following values in the following order:

        - Mean wind model reference speed.
        - Mean wind model reference direction.
    """
    ref_speed = get_tmp_name(suffix='.tif')
    ref_direction = get_tmp_name(suffix='.tif')
    
    out = []
    for src, dst in zip([wm_ref_speed, wm_ref_direction], [ref_speed, ref_direction]):
        gdalwarp(src=src, dst=dst,
                 outputBounds=bounds,
                 dstSRS=f'EPSG:{epsg}',
                 xRes=resolution, yRes=resolution,
                 resampleAlg='bilinear')
        with Raster(dst) as ras:
            arr = ras.array()
            out.append(round(np.nanmean(arr), 2))
        os.remove(dst)
    
    return tuple(out)


def copy_src_meta(ard_dir, src_ids):
    """
    Copies the original metadata of the source scenes to the ARD product
    directory.

    Parameters
    ----------
    ard_dir: str
        A path pointing to the current ARD product directory.
    src_ids: list[pyroSAR.drivers.ID]
        List of :class:`~pyroSAR.drivers.ID` objects of all source scenes that overlap with the current MGRS tile.

    Returns
    -------
    None
    """
    for src_id in src_ids:
        source_dir = os.path.join(ard_dir, 'source')
        pid = re.match(src_id.pattern, os.path.basename(src_id.file)).group('productIdentifier')
        
        if src_id.scene.endswith('.zip'):
            base = os.path.basename(src_id.file)
            target = os.path.join(source_dir, pid)
            if not os.path.isdir(target):
                with zipfile.ZipFile(src_id.scene, 'r') as zip_ref:
                    zip_ref.extract(member=base + '/manifest.safe', path=source_dir)
                    annotation_files = [f for f in zip_ref.namelist() if base + '/annotation' in f]
                    zip_ref.extractall(members=annotation_files, path=source_dir)
                os.rename(os.path.join(source_dir, base), target)
        else:
            pid_dir = os.path.join(source_dir, pid)
            os.makedirs(pid_dir, exist_ok=True)
            shutil.copy(src=os.path.join(src_id.scene, 'manifest.safe'),
                        dst=os.path.join(pid_dir, 'manifest.safe'))
            shutil.copytree(src=os.path.join(src_id.scene, 'annotation'),
                            dst=os.path.join(pid_dir, 'annotation'),
                            dirs_exist_ok=True)


def find_in_annotation(annotation_dict: dict[str, etree.ElementTree], pattern, single=False, out_type='str'):
    """
    Search for a pattern in all XML annotation files provided and return a dictionary of results.

    Parameters
    ----------
    annotation_dict: dict
        A dict of annotation files in the form: {'swath ID': `lxml.etree._Element` object}
    pattern: str
        The pattern to search for in each annotation file.
    single: bool
        If True, the results found in each annotation file are expected to be the same and therefore only a single
        value will be returned instead of a dict. If the results differ, an error is raised. Default is False.
    out_type: str
        Output type to convert the results to. Can be one of the following:

        - 'str' (default)
        - 'float'
        - 'int'

    Returns
    -------
    out: dict
        A dictionary of the results containing a list for each of the annotation files. E.g.,
        {'swath ID': list[str or float or int]}
    """
    out = {}
    for s, a in annotation_dict.items():
        swaths = [x.text for x in a.findall('.//swathProcParams/swath')]
        items = a.findall(pattern)
        
        parent = items[0].getparent().tag
        if parent in ['azimuthProcessing', 'rangeProcessing']:
            for i, val in enumerate(items):
                out[swaths[i]] = val.text
        else:
            out[s] = [x.text for x in items]
            if len(out[s]) == 1:
                out[s] = out[s][0]
    
    def _convert(obj, type):
        if isinstance(obj, list):
            return [_convert(x, type) for x in obj]
        elif isinstance(obj, str):
            if type == 'float':
                return float(obj)
            if type == 'int':
                return int(obj)
    
    if out_type != 'str':
        for k, v in list(out.items()):
            out[k] = _convert(v, out_type)
    
    err_msg = 'Search result for pattern "{}" expected to be the same in all annotation files.'
    if single:
        val = list(out.values())[0]
        for k in out:
            if out[k] != val:
                raise RuntimeError(err_msg.format(pattern))
        if out_type != 'str':
            return _convert(val, out_type)
        else:
            return val
    else:
        return out


def get_osv_info(sid):
    """
    Get information about the used OSV file.
    First, this function attempts to find an auxiliary OSV file matching the scene.
    If found, its name is returned. If not, it is assumed that it is not yet available
    and processing was performed using the OSVs found in the source product.
    In this case, the metadata is searched for the name of an auxiliary OSV file
    used during L1 generation. If found, its name is returned.

    Parameters
    ----------
    sid: pyroSAR.drivers.ID
        The pyroSAR scene ID object

    Returns
    -------
    tuple[str or None]
        the OSV file's basename and the OSV type description.
        None is returned if no OSV file is found.

    See Also
    --------
    pyroSAR.drivers.SAFE.getOSV
    """
    # try to find external OSV files
    osv = sid.getOSV(returnMatch=True, osvType=['POE', 'RES'], useLocal=True)
    if osv is None:
        # read the OSV file used during preprocessing from the metadata
        with sid.getFileObj(sid.findfiles('manifest.safe')[0]) as f:
            manifest = f.getvalue()
        tree = etree.fromstring(manifest)
        pattern = ("//safe:resource[@role='AUX_POE' or "
                   "@role='AUX_RES' or @role='AUX_PRE']")
        osv_match = tree.xpath(pattern, namespaces=tree.nsmap)
        if len(osv_match) > 0:
            osv = osv_match[0].get('name')
    osv_descr = None
    if osv is not None:
        while '.' in osv:
            osv = os.path.splitext(os.path.basename(osv))[0]
        osv_type = re.search('(?:POE|RES|PRE)ORB', osv).group()
        osv_descr = OSV_MAP[osv_type]
    return osv, osv_descr


def get_prod_meta(tif, src_ids, sar_dir, processor_name):
    """
    Returns a metadata dictionary, which is generated from the name of a product scene using a regular expression
    pattern and from a measurement GeoTIFF file of the same product scene using the :class:`~spatialist.raster.Raster`
    class.

    Parameters
    ----------
    tif: str
        The path to a measurement GeoTIFF file of the product scene.
    src_ids: list[pyroSAR.drivers.ID]
        List of :class:`~pyroSAR.drivers.ID` objects of all source SLC scenes that overlap with the current MGRS tile.
    sar_dir: str
        A path pointing to the processed SAR datasets of the product.
    processor_name: str
        The name of the SAR processor. Needed for reading processing metadata.

    Returns
    -------
    dict
        A dictionary containing metadata for the product scene.
    """
    processor = load_processor(processor_name)
    out = dict()
    coord_list = [sid.meta['coordinates'] for sid in src_ids]
    
    with Raster(tif) as ras:
        vec = ras.bbox()
        srs = vec.srs
        out['wkt'] = srs.ExportToWkt()
        out['epsg'] = vec.getProjection(type='epsg')
        out['rows'] = ras.rows
        out['cols'] = ras.cols
        out['res'] = ras.res
        geo = ras.geo
        out['transform'] = [geo['xres'], geo['rotation_x'], geo['xmin'],
                            geo['rotation_y'], geo['yres'], geo['ymax']]
        out['geom'] = geometry_from_vec(vectorobject=vec)
        
        # Calculate number of nodata border pixels based on source scene(s) footprint
        with vec_from_srccoords(coord_list=coord_list, crs=4326) as srcvec:
            ras_srcvec = rasterize(vectorobject=srcvec, reference=ras, burn_values=[1])
            arr_srcvec = ras_srcvec.array()
            out['nodata_borderpx'] = np.count_nonzero(np.isnan(arr_srcvec))
    
    proc_meta = processor.get_metadata(scene=src_ids[0].scene, outdir=sar_dir)
    out['ML_nRgLooks'] = proc_meta['rlks']
    out['ML_nAzLooks'] = proc_meta['azlks']
    return out


# ElementTree is a cython class, which may behave like a function and may
# thus not define __or__, which is needed for the | typing operator.
# Using Union instead.
def get_src_meta(
        sid: ID
) -> dict[str, Union[etree.ElementTree, dict[str, etree.ElementTree]]]:
    """
    Retrieve the manifest and annotation XML data of a scene as a dictionary
    using an :class:`pyroSAR.drivers.ID` object.

    Parameters
    ----------
    sid:
        A pyroSAR :class:`~pyroSAR.drivers.ID` object generated with e.g.
        :func:`pyroSAR.drivers.identify`.

    Returns
    -------
        A dictionary containing the parsed `etree.ElementTree` objects for
        the manifest and annotation XML files.
    """
    files = sid.findfiles(r'^s1[abcd].*-[vh]{2}-.*\.xml$')
    pols = list(set([re.search('[vh]{2}', os.path.basename(a)).group() for a in files]))
    annotation_files = list(filter(re.compile(pols[0]).search, files))
    
    a_files_base = [os.path.basename(a) for a in annotation_files]
    swaths = [re.search('-(iw[1-3]*|ew[1-5]*|s[1-6])', a).group(1) for a in a_files_base]
    
    annotation_dict = {}
    for s, a in zip(swaths, annotation_files):
        annotation_dict[s.upper()] = etree.fromstring(sid.getFileObj(a).getvalue())
    
    with sid.getFileObj(sid.findfiles('manifest.safe')[0]) as input_man:
        manifest = etree.fromstring(input_man.getvalue())
    
    return {'manifest': manifest,
            'annotation': annotation_dict}


def meta_dict(config, prod_meta, src_ids, compression):
    """
    Creates a dictionary containing metadata for a product scene, as well
    as its source scenes. The dictionary can then be used
    by :func:`~s1ard.metadata.xml.parse` and :func:`~s1ard.metadata.stac.parse`
    to generate OGC XML and STAC JSON metadata files, respectively.
    
    Parameters
    ----------
    config: dict
        Dictionary of the parsed config parameters for the current process.
    prod_meta: dict
        a metadata dictionary as returned by :func:`s1ard.ard.product_info`
    src_ids: list[pyroSAR.drivers.ID]
        List of :class:`~pyroSAR.drivers.ID` objects of all source scenes that overlap with the current MGRS tile.
    compression: str
        The compression type applied to raster files of the product.
    
    Returns
    -------
    meta: dict
        A dictionary containing a collection of metadata for product as well as source scenes.
    """
    dummy_num = -99999
    dummy_str = 'TBD'
    
    meta = {'prod': {},
            'source': {},
            'common': {}}
    src_sid = {}
    src_xml = {}
    for i, sid in enumerate(src_ids):
        uid = os.path.basename(sid.scene).split('.')[0][-4:]
        src_sid[uid] = sid
        src_xml[uid] = get_src_meta(sid=sid)
    sid0 = src_sid[list(src_sid.keys())[0]]  # first key/first file; used to extract some common metadata
    
    ref_tif = finder(prod_meta['dir_ard_product'], ['[hv]{2}-[gs]-lin.tif$'], regex=True)[0]
    np_tifs = finder(prod_meta['dir_ard_product'], ['-np-[hv]{2}.tif$'], regex=True)
    ei_tif = finder(prod_meta['dir_ard_product'], ['-ei.tif$'], regex=True)
    prod_meta.update(get_prod_meta(tif=ref_tif, src_ids=src_ids,
                                   sar_dir=config['processing']['sar_dir'],
                                   processor_name=config['processing']['processor']))
    op_mode = prod_meta['mode']
    
    # COMMON metadata (sorted alphabetically)
    meta['common']['antennaLookDirection'] = 'RIGHT'
    meta['common']['constellation'] = 'sentinel-1'
    meta['common']['instrumentShortName'] = 'C-SAR'
    meta['common']['operationalMode'] = op_mode
    meta['common']['orbitDirection'] = {'A': 'ascending', 'D': 'descending'}[sid0.orbit]
    meta['common']['orbitMeanAltitude'] = '{:.2e}'.format(693000)
    meta['common']['orbitNumber_abs'] = sid0.meta['orbitNumber_abs']
    meta['common']['orbitNumber_rel'] = sid0.meta['orbitNumber_rel']
    pid_lookup = {'S1A': 'A', 'S1B': 'B', 'S1C': 'C', 'S1D': 'D'}
    meta['common']['platformIdentifier'] = pid_lookup[sid0.sensor]
    meta['common']['platformShortName'] = 'Sentinel-1'
    meta['common']['platformFullname'] = '{}{}'.format(meta['common']['platformShortName'].lower(),
                                                       meta['common']['platformIdentifier'].lower())
    meta['common']['platformReference'] = URL['platformReference'][meta['common']['platformFullname']]
    meta['common']['polarisationChannels'] = sid0.polarizations
    meta['common']['polarisationMode'] = prod_meta['polarization'][0]
    meta['common']['processingLevel'] = 'L1C'
    meta['common']['radarBand'] = 'C'
    meta['common']['radarCenterFreq'] = 5405000000
    meta['common']['sensorType'] = 'RADAR'
    meta['common']['swathIdentifier'] = op_mode
    meta['common']['wrsLongitudeGrid'] = str(sid0.meta['orbitNumbers_rel']['start'])
    
    # Product metadata (sorted alphabetically)
    meta['prod']['access'] = config['metadata']['access_url']
    meta['prod']['acquisitionType'] = 'NOMINAL'
    meta['prod']['ancillaryData_KML'] = URL['ancillaryData_KML']
    meta['prod']['azimuthNumberOfLooks'] = round(prod_meta['ML_nAzLooks'], 2)
    meta['prod']['backscatterConvention'] = 'linear power'
    meta['prod']['backscatterConversionEq'] = '10*log10(DN)'
    meta['prod']['backscatterMeasurement'] = 'gamma0' if re.search('g-lin', ref_tif) else 'sigma0'
    if prod_meta['product_type'] == 'ORB':
        meta['prod']['card4l-link'] = URL['card4l_orb']
        meta['prod']['card4l-version'] = '1.0'
    else:
        meta['prod']['card4l-link'] = URL['card4l_nrb']
        meta['prod']['card4l-version'] = '5.5'
    meta['prod']['compression_type'] = compression
    meta['prod']['compression_zerrors'] = LERC_ERR_THRES
    meta['prod']['crsEPSG'] = str(prod_meta['epsg'])
    meta['prod']['crsWKT'] = prod_meta['wkt']
    meta['prod']['demAccess'] = DEM_MAP[config['processing']['dem_type']]['access']
    meta['prod']['demEGMReference'] = DEM_MAP[config['processing']['dem_type']]['egm']
    meta['prod']['demEGMResamplingMethod'] = 'bilinear'
    meta['prod']['demGSD'] = DEM_MAP[config['processing']['dem_type']]['gsd']
    meta['prod']['demName'] = config['processing']['dem_type'].replace(' II', '')
    meta['prod']['demReference'] = DEM_MAP[config['processing']['dem_type']]['ref']
    meta['prod']['demResamplingMethod'] = 'bilinear'
    meta['prod']['demType'] = DEM_MAP[config['processing']['dem_type']]['type']
    meta['prod']['doi'] = config['metadata']['doi']
    meta['prod']['ellipsoidalHeight'] = None
    meta['prod']['equivalentNumberOfLooks'] = calc_enl(tif=ref_tif)
    
    if (len(ei_tif) == 1 and
            sid0.product == 'SLC' and
            'copernicus' in config['processing']['dem_type'].lower()):
        geo_corr_accuracy = calc_geolocation_accuracy(swath_identifier=op_mode,
                                                      ei_tif=ei_tif[0],
                                                      etad=config['processing']['etad'])
    else:
        geo_corr_accuracy = None
    meta['prod']['geoCorrAccuracyEasternBias'] = dummy_num
    meta['prod']['geoCorrAccuracyEasternSTDev'] = dummy_num
    meta['prod']['geoCorrAccuracyNorthernBias'] = dummy_num
    meta['prod']['geoCorrAccuracyNorthernSTDev'] = dummy_num
    if geo_corr_accuracy is not None:
        meta['prod']['geoCorrAccuracyReference'] = URL['geoCorrAccuracyReference']
        meta['prod']['geoCorrAccuracy_rRMSE'] = geo_corr_accuracy
    else:
        meta['prod']['geoCorrAccuracyReference'] = dummy_str
        meta['prod']['geoCorrAccuracy_rRMSE'] = dummy_num
    meta['prod']['geoCorrAccuracyType'] = 'gtc'  # or slant-range
    
    meta['prod']['geoCorrAlgorithm'] = URL['geoCorrAlgorithm']
    meta['prod']['geoCorrResamplingMethod'] = 'bilinear'
    meta['prod']['geom_stac_bbox_native'] = prod_meta['geom']['bbox_native']
    meta['prod']['geom_stac_bbox_4326'] = prod_meta['geom']['bbox']
    meta['prod']['geom_stac_geometry_4326'] = prod_meta['geom']['geometry']
    meta['prod']['geom_xml_center'] = prod_meta['geom']['center']
    meta['prod']['geom_xml_envelope'] = prod_meta['geom']['envelope']
    meta['prod']['griddingConvention'] = 'Military Grid Reference System (MGRS)'
    meta['prod']['griddingConventionURL'] = URL['griddingConventionURL']
    meta['prod']['licence'] = config['metadata']['licence']
    meta['prod']['mgrsID'] = prod_meta['tile']
    meta['prod']['noiseRemovalApplied'] = True
    nr_algo = URL['noiseRemovalAlgorithm'] if meta['prod']['noiseRemovalApplied'] else None
    meta['prod']['noiseRemovalAlgorithm'] = nr_algo
    meta['prod']['numberOfAcquisitions'] = str(len(src_sid))
    meta['prod']['numBorderPixels'] = prod_meta['nodata_borderpx']
    meta['prod']['numLines'] = str(prod_meta['rows'])
    meta['prod']['numPixelsPerLine'] = str(prod_meta['cols'])
    meta['prod']['pixelCoordinateConvention'] = 'upper-left'
    processing_center = config['metadata']['processing_center']
    if processing_center is None:
        processing_center = 'None'
    meta['prod']['processingCenter'] = processing_center
    meta['prod']['processingMode'] = 'PROTOTYPE'
    meta['prod']['processorName'] = 's1ard'
    meta['prod']['processorVersion'] = s1ard.__version__
    prod_name_prefix = 'Ocean' if prod_meta['product_type'] == 'ORB' else 'Normalised'
    meta['prod']['productName'] = f"{prod_name_prefix} Radar Backscatter"
    meta['prod']['productName-short'] = prod_meta['product_type']
    meta['prod']['pxSpacingColumn'] = str(prod_meta['res'][0])
    meta['prod']['pxSpacingRow'] = str(prod_meta['res'][1])
    meta['prod']['radiometricAccuracyAbsolute'] = dummy_num
    meta['prod']['radiometricAccuracyRelative'] = dummy_num
    meta['prod']['radiometricAccuracyReference'] = URL['radiometricAccuracyReference']
    meta['prod']['rangeNumberOfLooks'] = round(prod_meta['ML_nRgLooks'], 2)
    meta['prod']['RTCAlgorithm'] = URL['RTCAlgorithm']
    meta['prod']['speckleFilterApplied'] = None
    meta['prod']['status'] = 'PLANNED'
    meta['prod']['timeCreated'] = prod_meta['proc_time']
    meta['prod']['timeStart'] = prod_meta['start']
    meta['prod']['timeStop'] = prod_meta['stop']
    meta['prod']['transform'] = prod_meta['transform']
    
    # SOURCE metadata
    for uid in list(src_sid.keys()):
        sid = src_sid[uid]
        nsmap = src_xml[uid]['manifest'].nsmap
        
        swath_ids = find_in_annotation(annotation_dict=src_xml[uid]['annotation'],
                                       pattern='.//swathProcParams/swath')
        swaths = []
        for item in swath_ids.values():
            if isinstance(item, list):
                swaths.extend(item)
            else:
                swaths.append(item)
        with sid.geometry() as vec:
            geom = geometry_from_vec(vectorobject=vec)
        
        res_mode = re.match(re.compile(sid.pattern), os.path.basename(sid.file)).groupdict()['resolution']
        product_type = sid.meta['product'] + (res_mode if res_mode != '_' else '')
        if product_type.startswith('GRD'):
            data_geometry = 'ground-range'
            if re.search('S[1-6]', op_mode):
                res_az = {op_mode: RES_MAP_GRD[res_mode]['SM']['az'][op_mode]}
                res_rg = {op_mode: RES_MAP_GRD[res_mode]['SM']['rg'][op_mode]}
                enl = round(mean(ENL_MAP_GRD[res_mode]['SM']), 2)
            else:
                res_az = RES_MAP_GRD[res_mode][op_mode]['az']
                res_rg = RES_MAP_GRD[res_mode][op_mode]['rg']
                enl = round(mean(ENL_MAP_GRD[res_mode][op_mode]), 2)
        else:  # SLC
            data_geometry = 'slant-range'
            res_az = RES_MAP_SLC[op_mode]['az']
            res_rg = RES_MAP_SLC[op_mode]['rg']
            enl = 1.0
        
        patterns = [
            './/azimuthProcessing/lookBandwidth',
            './/rangeProcessing/lookBandwidth',
            './/azimuthProcessing/numberOfLooks',
            './/rangeProcessing/numberOfLooks',
            './/azimuthPixelSpacing',
            './/rangePixelSpacing',
            './/geolocationGridPoint/incidenceAngle'
        ]
        out_types = ['float', 'float', 'int', 'int', 'float', 'float', 'float']
        results = []
        for pattern, out_type in zip(patterns, out_types):
            result = find_in_annotation(annotation_dict=src_xml[uid]['annotation'],
                                        pattern=pattern, out_type=out_type)
            results.append(result)
        az_look_bandwidth, rg_look_bandwidth, az_num_looks, rg_num_looks, az_px_spacing, rg_px_spacing, inc = results
        
        inc_vals = dissolve(list(inc.values()))
        lut_applied = find_in_annotation(annotation_dict=src_xml[uid]['annotation'],
                                         pattern='.//applicationLutId', single=True)
        pslr, islr = calc_pslr_islr(annotation_dict=src_xml[uid]['annotation'])
        
        def _read_manifest(pattern, attrib=None):
            obj = src_xml[uid]['manifest'].find(pattern, nsmap)
            if attrib is not None:
                return obj.attrib[attrib]
            else:
                return obj.text
        
        # Source product metadata (sorted alphabetically)
        meta['source'][uid] = {}
        meta['source'][uid]['access'] = URL['source_access']
        meta['source'][uid]['acquisitionType'] = 'NOMINAL'
        asc_node_time = dateparse(_read_manifest('.//s1:ascendingNodeTime'))
        asc_node_time = asc_node_time.replace(tzinfo=timezone.utc)
        meta['source'][uid]['ascendingNodeDate'] = asc_node_time
        meta['source'][uid]['azimuthLookBandwidth'] = az_look_bandwidth
        meta['source'][uid]['azimuthNumberOfLooks'] = az_num_looks
        meta['source'][uid]['azimuthPixelSpacing'] = az_px_spacing
        meta['source'][uid]['azimuthResolution'] = res_az
        meta['source'][uid]['dataGeometry'] = data_geometry
        meta['source'][uid]['datatakeID'] = _read_manifest('.//s1sarl1:missionDataTakeID')
        meta['source'][uid]['doi'] = URL['source_doi']
        meta['source'][uid]['faradayMeanRotationAngle'] = None
        meta['source'][uid]['faradayRotationReference'] = URL['faradayRotationReference']
        meta['source'][uid]['filename'] = sid.file
        meta['source'][uid]['geom_stac_bbox_4326'] = geom['bbox']
        meta['source'][uid]['geom_stac_geometry_4326'] = geom['geometry']
        meta['source'][uid]['geom_xml_center'] = geom['center']
        meta['source'][uid]['geom_xml_envelope'] = geom['envelope']
        meta['source'][uid]['incidenceAngleMax'] = round(np.max(inc_vals), 2)
        meta['source'][uid]['incidenceAngleMin'] = round(np.min(inc_vals), 2)
        meta['source'][uid]['incidenceAngleMidSwath'] = round(np.max(inc_vals) -
                                                              ((np.max(inc_vals) - np.min(inc_vals)) / 2), 2)
        meta['source'][uid]['instrumentAzimuthAngle'] = round(sid.meta['heading'], 2)
        meta['source'][uid]['ionosphereIndicator'] = None
        meta['source'][uid]['lutApplied'] = lut_applied
        meta['source'][uid]['majorCycleID'] = str(sid.meta['cycleNumber'])
        meta['source'][uid]['orbitDataAccess'] = URL['orbitDataAccess']
        osv_base, osv_descr = get_osv_info(sid)
        meta['source'][uid]['orbitStateVector'] = osv_base
        meta['source'][uid]['orbitDataSource'] = osv_descr
        if len(np_tifs) > 0:
            meta['source'][uid]['perfEstimates'] = calc_performance_estimates(files=np_tifs)
            meta['source'][uid]['perfNoiseEquivalentIntensityType'] = 'sigma0'
        else:
            stats = {stat: None for stat in ['minimum', 'mean', 'maximum']}
            pe = {pol: stats for pol in meta['common']['polarisationChannels']}
            meta['source'][uid]['perfEstimates'] = pe
            meta['source'][uid]['perfNoiseEquivalentIntensityType'] = None
        meta['source'][uid]['perfEquivalentNumberOfLooks'] = enl
        meta['source'][uid]['perfIntegratedSideLobeRatio'] = islr
        meta['source'][uid]['perfPeakSideLobeRatio'] = pslr
        meta['source'][uid]['polCalMatrices'] = None
        fac_org = _read_manifest('.//safe:facility', attrib='organisation')
        fac_name = _read_manifest('.//safe:facility', attrib='name')
        meta['source'][uid]['processingCenter'] = f"{fac_org} {fac_name}".replace(' -', '')
        proc_date = dateparse(_read_manifest('.//safe:processing', attrib='stop'))
        proc_date = proc_date.replace(tzinfo=timezone.utc)
        meta['source'][uid]['processingDate'] = proc_date
        meta['source'][uid]['processingLevel'] = _read_manifest('.//safe:processing', attrib='name')
        meta['source'][uid]['processingMode'] = 'NOMINAL'
        meta['source'][uid]['processorName'] = _read_manifest('.//safe:software', attrib='name')
        meta['source'][uid]['processorVersion'] = _read_manifest('.//safe:software', attrib='version')
        meta['source'][uid]['productType'] = product_type
        meta['source'][uid]['rangeLookBandwidth'] = rg_look_bandwidth
        meta['source'][uid]['rangeNumberOfLooks'] = rg_num_looks
        meta['source'][uid]['rangePixelSpacing'] = rg_px_spacing
        meta['source'][uid]['rangeResolution'] = res_rg
        meta['source'][uid]['sensorCalibration'] = URL['sensorCalibration']
        meta['source'][uid]['status'] = 'ARCHIVED'
        meta['source'][uid]['swaths'] = swaths
        meta['source'][uid]['timeCompletionFromAscendingNode'] = str(float(_read_manifest('.//s1:stopTimeANX')))
        meta['source'][uid]['timeStartFromAscendingNode'] = str(float(_read_manifest('.//s1:startTimeANX')))
        meta['source'][uid]['timeStart'] = dateparse(sid.start).replace(tzinfo=timezone.utc)
        meta['source'][uid]['timeStop'] = dateparse(sid.stop).replace(tzinfo=timezone.utc)
    return meta
