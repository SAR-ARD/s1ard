import os
import re
import shutil
import zipfile
import math
from statistics import mean
import json
from lxml import etree
from dateutil.parser import parse as dateparse
import numpy as np
from statistics import median
from spatialist import Raster
from spatialist.auxil import gdalwarp, crsConvert
from spatialist.ancillary import finder, dissolve
from spatialist.raster import rasterize
from spatialist.vector import Vector
from osgeo import gdal, ogr
import s1ard
from s1ard.metadata.mapping import (ARD_PATTERN, LERC_ERR_THRES, RES_MAP_SLC, RES_MAP_GRD,
                                    ENL_MAP_GRD, OSV_MAP, DEM_MAP, SLC_ACC_MAP, URL)
from s1ard import snap
from s1ard.ancillary import get_tmp_name

gdal.UseExceptions()


def meta_dict(config, target, src_ids, sar_dir, proc_time, start, stop,
              compression, product_type, wm_ref_files=None):
    """
    Creates a dictionary containing metadata for a product scene, as well
    as its source scenes. The dictionary can then be used
    by :func:`~s1ard.metadata.xml.parse` and :func:`~s1ard.metadata.stac.parse`
    to generate OGC XML and STAC JSON metadata files, respectively.
    
    Parameters
    ----------
    config: dict
        Dictionary of the parsed config parameters for the current process.
    target: str
        A path pointing to the current ARD product directory.
    src_ids: list[pyroSAR.drivers.ID]
        List of :class:`~pyroSAR.drivers.ID` objects of all source scenes that overlap with the current MGRS tile.
    sar_dir: str
        The SAR processing output directory.
    proc_time: datetime.datetime
        The processing time object used to generate the unique product identifier.
    start: datetime.datetime
        The product start time.
    stop: datetime.datetime
        The product stop time.
    compression: str
        The compression type applied to raster files of the product.
    product_type: str
        The type of ARD product that is being created. Either 'NRB' or 'ORB'.
    wm_ref_files: list[str], optional
        A list of paths pointing to wind model reference files. Default is None.
    
    Returns
    -------
    meta: dict
        A dictionary containing a collection of metadata for product as well as source scenes.
    """
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
    swath_id = re.search('_(IW|EW|S[1-6])_', os.path.basename(sid0.file)).group().replace('_', '')
    
    ref_tif = finder(target, ['[hv]{2}-[gs]-lin.tif$'], regex=True)[0]
    ratio_tif = finder(target, ['[hv]{2}-[gs]-lin.vrt$'], regex=True)
    np_tifs = finder(target, ['-np-[hv]{2}.tif$'], regex=True)
    ei_tif = finder(target, ['-ei.tif$'], regex=True)
    product_id = os.path.basename(target)
    prod_meta = get_prod_meta(product_id=product_id, tif=ref_tif,
                              src_ids=src_ids, sar_dir=sar_dir)
    op_mode = prod_meta['mode']
    
    tups = [(key, LERC_ERR_THRES[key]) for key in LERC_ERR_THRES.keys()]
    z_err_dict = dict(tups)
    
    # COMMON metadata (sorted alphabetically)
    meta['common']['antennaLookDirection'] = 'RIGHT'
    meta['common']['constellation'] = 'sentinel-1'
    meta['common']['instrumentShortName'] = 'C-SAR'
    meta['common']['operationalMode'] = op_mode
    meta['common']['orbitDirection'] = {'A': 'ascending', 'D': 'descending'}[sid0.orbit]
    meta['common']['orbitMeanAltitude'] = '{:.2e}'.format(693000)
    meta['common']['orbitNumber_abs'] = sid0.meta['orbitNumber_abs']
    meta['common']['orbitNumber_rel'] = sid0.meta['orbitNumber_rel']
    pid_lookup = {'S1A': '1A', 'S1B': '1B', 'S1C': '1C', 'S1D': '1D'}
    meta['common']['platformIdentifier'] = pid_lookup[sid0.sensor]
    meta['common']['platformShortName'] = 'Sentinel'
    meta['common']['platformFullname'] = '{}-{}'.format(meta['common']['platformShortName'].lower(),
                                                        meta['common']['platformIdentifier'].lower())
    meta['common']['platformReference'] = URL['platformReference'][meta['common']['platformFullname']]
    meta['common']['polarisationChannels'] = sid0.polarizations
    meta['common']['polarisationMode'] = prod_meta['pols'][0]
    meta['common']['processingLevel'] = 'L1C'
    meta['common']['radarBand'] = 'C'
    meta['common']['radarCenterFreq'] = 5405000000
    meta['common']['sensorType'] = 'RADAR'
    meta['common']['swathIdentifier'] = swath_id
    meta['common']['wrsLongitudeGrid'] = str(sid0.meta['orbitNumbers_rel']['start'])
    
    # PRODUCT metadata
    if (len(ei_tif) == 1 and
            sid0.product == 'SLC' and
            'copernicus' in config['processing']['dem_type'].lower()):
        geo_corr_accuracy = calc_geolocation_accuracy(swath_identifier=swath_id,
                                                      ei_tif=ei_tif[0],
                                                      etad=config['processing']['etad'])
    else:
        geo_corr_accuracy = None
    
    # (sorted alphabetically)
    meta['prod']['access'] = config['metadata']['access_url']
    meta['prod']['acquisitionType'] = 'NOMINAL'
    meta['prod']['ancillaryData_KML'] = URL['ancillaryData_KML']
    meta['prod']['azimuthNumberOfLooks'] = round(prod_meta['ML_nAzLooks'], 2)
    meta['prod']['backscatterConvention'] = 'linear power'
    meta['prod']['backscatterConversionEq'] = '10*log10(DN)'
    meta['prod']['backscatterMeasurement'] = 'gamma0' if re.search('g-lin', ref_tif) else 'sigma0'
    if product_type == 'ORB':
        meta['prod']['card4l-link'] = URL['card4l_orb']
        meta['prod']['card4l-version'] = '1.0'
    else:
        meta['prod']['card4l-link'] = URL['card4l_nrb']
        meta['prod']['card4l-version'] = '5.5'
    meta['prod']['compression_type'] = compression
    meta['prod']['compression_zerrors'] = z_err_dict
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
    meta['prod']['equivalentNumberLooks'] = calc_enl(tif=ref_tif)
    meta['prod']['geoCorrAccuracyEasternBias'] = None
    meta['prod']['geoCorrAccuracyEasternSTDev'] = None
    meta['prod']['geoCorrAccuracyNorthernBias'] = None
    meta['prod']['geoCorrAccuracyNorthernSTDev'] = None
    if geo_corr_accuracy is not None:
        geo_corr_accuracy_reference = URL['geoCorrAccuracyReference']
        geo_corr_accuracy_type = 'slant-range'
    else:
        geo_corr_accuracy_reference = None
        geo_corr_accuracy_type = None
    meta['prod']['geoCorrAccuracyReference'] = geo_corr_accuracy_reference
    meta['prod']['geoCorrAccuracyType'] = geo_corr_accuracy_type
    meta['prod']['geoCorrAccuracy_rRMSE'] = geo_corr_accuracy
    
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
    meta['prod']['mgrsID'] = prod_meta['mgrsID']
    meta['prod']['noiseRemovalApplied'] = True
    nr_algo = URL['noiseRemovalAlgorithm'] if meta['prod']['noiseRemovalApplied'] else None
    meta['prod']['noiseRemovalAlgorithm'] = nr_algo
    meta['prod']['numberOfAcquisitions'] = str(len(src_sid))
    meta['prod']['numBorderPixels'] = prod_meta['nodata_borderpx']
    meta['prod']['numLines'] = str(prod_meta['rows'])
    meta['prod']['numPixelsPerLine'] = str(prod_meta['cols'])
    meta['prod']['pixelCoordinateConvention'] = 'upper-left'
    meta['prod']['processingCenter'] = config['metadata']['processing_center']
    meta['prod']['processingMode'] = 'PROTOTYPE'
    meta['prod']['processorName'] = 's1ard'
    meta['prod']['processorVersion'] = s1ard.__version__
    meta['prod']['productName'] = f"{'Ocean' if product_type == 'ORB' else 'Normalised'} Radar Backscatter"
    meta['prod']['productName-short'] = product_type
    meta['prod']['pxSpacingColumn'] = str(prod_meta['res'][0])
    meta['prod']['pxSpacingRow'] = str(prod_meta['res'][1])
    meta['prod']['radiometricAccuracyAbsolute'] = None
    meta['prod']['radiometricAccuracyRelative'] = None
    meta['prod']['radiometricAccuracyReference'] = URL['radiometricAccuracyReference']
    meta['prod']['rangeNumberOfLooks'] = round(prod_meta['ML_nRgLooks'], 2)
    meta['prod']['RTCAlgorithm'] = URL['RTCAlgorithm']
    meta['prod']['speckleFilterApplied'] = False
    meta['prod']['status'] = 'PLANNED'
    meta['prod']['timeCreated'] = proc_time
    meta['prod']['timeStart'] = start
    meta['prod']['timeStop'] = stop
    meta['prod']['transform'] = prod_meta['transform']
    if wm_ref_files is not None:
        wm_ref_mean_speed, wm_ref_mean_dir = calc_wm_ref_stats(wm_ref_files=wm_ref_files,
                                                               epsg=prod_meta['epsg'],
                                                               bounds=prod_meta['geom']['bbox_native'])
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
        meta['source'][uid]['ascendingNodeDate'] = _read_manifest('.//s1:ascendingNodeTime')
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
        meta['source'][uid]['geom_xml_envelop'] = geom['envelop']
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
        meta['source'][uid]['timeStart'] = dateparse(sid.start)
        meta['source'][uid]['timeStop'] = dateparse(sid.stop)
    
    return meta


def get_prod_meta(product_id, tif, src_ids, sar_dir):
    """
    Returns a metadata dictionary, which is generated from the name of a product scene using a regular expression
    pattern and from a measurement GeoTIFF file of the same product scene using the :class:`~spatialist.raster.Raster`
    class.
    
    Parameters
    ----------
    product_id: str
        The top-level product folder name.
    tif: str
        The path to a measurement GeoTIFF file of the product scene.
    src_ids: list[pyroSAR.drivers.ID]
        List of :class:`~pyroSAR.drivers.ID` objects of all source SLC scenes that overlap with the current MGRS tile.
    sar_dir: str
        A path pointing to the processed SAR datasets of the product.
    
    Returns
    -------
    dict
        A dictionary containing metadata for the product scene.
    """
    out = re.match(re.compile(ARD_PATTERN), product_id).groupdict()
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
        with _vec_from_srccoords(coord_list=coord_list, crs=4326) as srcvec:
            ras_srcvec = rasterize(vectorobject=srcvec, reference=ras, burn_values=[1])
            arr_srcvec = ras_srcvec.array()
            out['nodata_borderpx'] = np.count_nonzero(np.isnan(arr_srcvec))
    
    src_xml = get_src_meta(sid=src_ids[0])
    az_num_looks = find_in_annotation(annotation_dict=src_xml['annotation'],
                                      pattern='.//azimuthProcessing/numberOfLooks',
                                      out_type='int')
    rg_num_looks = find_in_annotation(annotation_dict=src_xml['annotation'],
                                      pattern='.//rangeProcessing/numberOfLooks',
                                      out_type='int')
    proc_meta = snap.get_metadata(scene=src_ids[0].scene, outdir=sar_dir)
    out['ML_nRgLooks'] = proc_meta['rlks'] * median(rg_num_looks.values())
    out['ML_nAzLooks'] = proc_meta['azlks'] * median(az_num_looks.values())
    return out


def _vec_from_srccoords(coord_list, crs, layername='polygon'):
    """
    Creates a single :class:`~spatialist.vector.Vector` object from a list
    of footprint coordinates of source scenes.
    
    Parameters
    ----------
    coord_list: list[list[tuple[float]]]
        List containing for each source scene a list of coordinate pairs as
        retrieved from the metadata stored in an :class:`~pyroSAR.drivers.ID`
        object.
    crs: int or str
        the coordinate reference system of the provided coordinates.
    layername: str
        the layer name of the output vector object
    
    Returns
    -------
    spatialist.vector.Vector
    """
    srs = crsConvert(crs, 'osr')
    pts = ogr.Geometry(ogr.wkbMultiPoint)
    for footprint in coord_list:
        for lon, lat in footprint:
            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(lon, lat)
            pts.AddGeometry(point)
    geom = pts.ConvexHull()
    vec = Vector(driver='Memory')
    vec.addlayer(layername, srs, geom.GetGeometryType())
    vec.addfeature(geom)
    point = None
    pts = None
    geom = None
    return vec


def get_src_meta(sid):
    """
    Retrieve the manifest and annotation XML data of a scene as a dictionary using an :class:`pyroSAR.drivers.ID`
    object.
    
    Parameters
    ----------
    sid:  pyroSAR.drivers.ID
        A pyroSAR :class:`~pyroSAR.drivers.ID` object generated with e.g. :func:`pyroSAR.drivers.identify`.
    
    Returns
    -------
    dict
        A dictionary containing the parsed `etree.ElementTree` objects for the manifest and annotation XML files.
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


def geometry_from_vec(vectorobject):
    """
    Get geometry information for usage in STAC and XML metadata from a :class:`spatialist.vector.Vector` object.
    
    Parameters
    ----------
    vectorobject: spatialist.vector.Vector
        The vector object to extract geometry information from.
    
    Returns
    -------
    out: dict
        A dictionary containing the geometry information extracted from the vector object.
    """
    out = {}
    vec = vectorobject
    
    # For STAC metadata
    if vec.getProjection(type='epsg') != 4326:
        ext = vec.extent
        out['bbox_native'] = [ext['xmin'], ext['ymin'], ext['xmax'], ext['ymax']]
        vec.reproject(4326)
    feat = vec.getfeatures()[0]
    geom = feat.GetGeometryRef()
    out['geometry'] = json.loads(geom.ExportToJson())
    ext = vec.extent
    out['bbox'] = [ext['xmin'], ext['ymin'], ext['xmax'], ext['ymax']]
    
    # For XML metadata
    c_x = (ext['xmax'] + ext['xmin']) / 2
    c_y = (ext['ymax'] + ext['ymin']) / 2
    out['center'] = '{} {}'.format(c_y, c_x)
    wkt = geom.ExportToWkt().removeprefix('POLYGON ((').removesuffix('))')
    wkt_list = ['{} {}'.format(x[1], x[0]) for x in [y.split(' ') for y in wkt.split(',')]]
    out['envelope'] = ' '.join(wkt_list)
    
    return out


def find_in_annotation(annotation_dict, pattern, single=False, out_type='str'):
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


def calc_enl(tif, block_size=30, return_arr=False, decimals=2):
    """
    Calculate the Equivalent Number of Looks (ENL) for a linear-scaled backscatter
    measurement GeoTIFF file. The calculation is performed block-wise for the
    entire image and by default the median ENL value is returned.
    
    Parameters
    ----------
    tif: str
        The path to a linear-scaled backscatter measurement GeoTIFF file.
    block_size: int, optional
        The block size to use for the calculation. Remainder pixels are discarded,
         if the array dimensions are not evenly divisible by the block size.
         Default is 30, which calculates ENL for 30x30 pixel blocks.
    return_arr: bool, optional
        If True, the calculated ENL array is returned. Default is False.
    decimals: int, optional
        Number of decimal places to round the calculated ENL value to. Default is 2.
    
    Raises
    ------
    RuntimeError
        if the input array contains only NaN values
    
    Returns
    -------
    float or None or numpy.ndarray
        If `return_enl_arr=True`, an array of ENL values is returned. Otherwise,
        the median ENL value is returned. If the ENL array contains only NaN and
        `return_enl_arr=False`, the return value is `None`.
    
    References
    ----------
    :cite:`anfinsen.etal_2009`
    """
    with Raster(tif) as ras:
        arr = ras.array()
    arr[np.isinf(arr)] = np.nan
    
    if len(arr[~np.isnan(arr)]) == 0:
        raise RuntimeError('cannot compute ENL for an empty array')
    
    num_blocks_rows = arr.shape[0] // block_size
    num_blocks_cols = arr.shape[1] // block_size
    if num_blocks_rows == 0 or num_blocks_cols == 0:
        raise ValueError("Block size is too large for the input data dimensions.")
    blocks = arr[:num_blocks_rows * block_size,
             :num_blocks_cols * block_size].reshape(num_blocks_rows, block_size,
                                                    num_blocks_cols, block_size)
    
    with np.testing.suppress_warnings() as sup:
        sup.filter(RuntimeWarning, "Mean of empty slice")
        _mean = np.nanmean(blocks, axis=(1, 3))
    with np.testing.suppress_warnings() as sup:
        sup.filter(RuntimeWarning, "Degrees of freedom <= 0 for slice")
        _std = np.nanstd(blocks, axis=(1, 3))
    enl = np.divide(_mean ** 2, _std ** 2,
                    out=np.full_like(_mean, fill_value=np.nan), where=_std != 0)
    out_arr = np.zeros((num_blocks_rows, num_blocks_cols))
    out_arr[:num_blocks_rows, :num_blocks_cols] = enl
    if not return_arr:
        if len(enl[~np.isnan(enl)]) == 0:
            return None
        out_median = np.nanmedian(out_arr)
        return np.round(out_median, decimals)
    else:
        return out_arr


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


def calc_performance_estimates(files, decimals=2):
    """
    Calculates the performance estimates specified in CARD4L NRB 1.6.9 for all noise power images if available.
    
    Parameters
    ----------
    files: list[str]
        List of paths pointing to the noise power images the estimates should be calculated for.
    decimals: int, optional
        Number of decimal places to round the calculated values to. Default is 2.
    
    Returns
    -------
    out: dict
        Dictionary containing the calculated estimates for each available polarization.
    """
    out = {}
    for f in files:
        pol = re.search('np-([vh]{2})', f).group(1).upper()
        with Raster(f) as ras:
            arr = ras.array()
            # The following need to be of type float, not numpy.float32 in order to be JSON serializable
            _min = float(np.nanmin(arr))
            _max = float(np.nanmax(arr))
            _mean = float(np.nanmean(arr))
            del arr
        out[pol] = {'minimum': round(_min, decimals),
                    'maximum': round(_max, decimals),
                    'mean': round(_mean, decimals)}
    return out


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


def calc_wm_ref_stats(wm_ref_files, epsg, bounds, resolution=915):
    """
    Calculates the mean wind model reference speed and direction for the wind model annotation layer.
    
    Parameters
    ----------
    wm_ref_files: list[str]
        List of paths pointing to the wind model reference files.
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
    files_speed = [f for f in wm_ref_files if f.endswith('Speed.tif')]
    files_direction = [f for f in wm_ref_files if f.endswith('Direction.tif')]
    
    ref_speed = get_tmp_name(suffix='.tif')
    ref_direction = get_tmp_name(suffix='.tif')
    
    out = []
    for src, dst in zip([files_speed, files_direction], [ref_speed, ref_direction]):
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


def get_header_size(tif):
    """
    Gets the header size of a GeoTIFF file in bytes.
    The code used in this function and its helper function `_get_block_offset` were extracted from the following
    source:
    
    https://github.com/OSGeo/gdal/blob/master/swig/python/gdal-utils/osgeo_utils/samples/validate_cloud_optimized_geotiff.py
    
    Copyright (c) 2017, Even Rouault
    
    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included
    in all copies or substantial portions of the Software.
    
    Parameters
    ----------
    tif: str
        A path to a GeoTIFF file of the currently processed ARD product.

    Returns
    -------
    header_size: int
        The size of all IFD headers of the GeoTIFF file in bytes.
    """
    
    def _get_block_offset(band):
        blockxsize, blockysize = band.GetBlockSize()
        for y in range(int((band.YSize + blockysize - 1) / blockysize)):
            for x in range(int((band.XSize + blockxsize - 1) / blockxsize)):
                block_offset = band.GetMetadataItem('BLOCK_OFFSET_%d_%d' % (x, y), 'TIFF')
                if block_offset:
                    return int(block_offset)
        return 0
    
    details = {}
    ds = gdal.Open(tif)
    main_band = ds.GetRasterBand(1)
    ovr_count = main_band.GetOverviewCount()
    
    block_offset = _get_block_offset(band=main_band)
    details['data_offsets'] = {}
    details['data_offsets']['main'] = block_offset
    for i in range(ovr_count):
        ovr_band = ds.GetRasterBand(1).GetOverview(i)
        block_offset = _get_block_offset(band=ovr_band)
        details['data_offsets']['overview_%d' % i] = block_offset
    
    headers_size = min(details['data_offsets'][k] for k in details['data_offsets'])
    if headers_size == 0:
        headers_size = gdal.VSIStatL(tif).size
    return headers_size


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
