import os
import re
import shutil
import zipfile
import math
from statistics import mean
import json
from lxml import etree
from datetime import datetime
import numpy as np
from statistics import median
from spatialist import Raster
from spatialist.ancillary import finder, dissolve
from spatialist.vector import wkt2vector
from spatialist.raster import rasterize
from osgeo import gdal
import S1_NRB
from S1_NRB.metadata.mapping import (ARD_PATTERN, LERC_ERR_THRES, RES_MAP_SLC, RES_MAP_GRD, ENL_MAP_GRD, OSV_MAP,
                                     DEM_MAP, SLC_ACC_MAP)
from S1_NRB import snap

gdal.UseExceptions()


def meta_dict(config, target, src_ids, sar_dir, proc_time, start, stop, compression, product_type):
    """
    Creates a dictionary containing metadata for a product scene, as well as its source scenes. The dictionary can then
    be utilized by :func:`~S1_NRB.metadata.xml.parse` and :func:`~S1_NRB.metadata.stac.parse` to generate OGC XML and
    STAC JSON metadata files, respectively.
    
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
    meta['common']['platformIdentifier'] = {'S1A': '1A', 'S1B': '1B'}[sid0.sensor]
    meta['common']['platformShortName'] = 'Sentinel'
    meta['common']['platformFullname'] = '{}-{}'.format(meta['common']['platformShortName'].lower(),
                                                        meta['common']['platformIdentifier'].lower())
    meta['common']['platformReference'] = \
        {'sentinel-1a': 'http://database.eohandbook.com/database/missionsummary.aspx?missionID=575',
         'sentinel-1b': 'http://database.eohandbook.com/database/missionsummary.aspx?missionID=576'}[
            meta['common']['platformFullname']]
    meta['common']['polarisationChannels'] = sid0.polarizations
    meta['common']['polarisationMode'] = prod_meta['pols'][0]
    meta['common']['processingLevel'] = 'L1C'
    meta['common']['radarBand'] = 'C'
    meta['common']['radarCenterFreq'] = 5405000000
    meta['common']['sensorType'] = 'RADAR'
    meta['common']['swathIdentifier'] = swath_id
    meta['common']['wrsLongitudeGrid'] = str(sid0.meta['orbitNumbers_rel']['start'])
    
    # PRODUCT metadata
    if len(ei_tif) == 1 and sid0.product == 'SLC' and 'copernicus' in config['dem_type'].lower():
        geo_corr_accuracy = calc_geolocation_accuracy(swath_identifier=swath_id, ei_tif=ei_tif[0], etad=config['etad'])
    else:
        geo_corr_accuracy = None
    
    # (sorted alphabetically)
    meta['prod']['access'] = config['meta']['access_url']
    meta['prod']['acquisitionType'] = 'NOMINAL'
    meta['prod']['ancillaryData_KML'] = 'https://sentinel.esa.int/documents/247904/1955685/S2A_OPER_GIP_TILPAR_MPC__' \
                                        '20151209T095117_V20150622T000000_21000101T000000_B00.kml'
    meta['prod']['azimuthNumberOfLooks'] = round(prod_meta['ML_nAzLooks'], 2)
    meta['prod']['backscatterConvention'] = 'linear power'
    meta['prod']['backscatterConversionEq'] = '10*log10(DN)'
    meta['prod']['backscatterMeasurement'] = 'gamma0' if re.search('g-lin', ref_tif) else 'sigma0'
    if product_type == 'ORB':
        meta['prod']['card4l-link'] = 'https://ceos.org/ard/files/PFS/ORB/v1.0/CARD4L_Product_Family_Specification_' \
                                      'Ocean_Radar_Backscatter-v1.0.pdf'
        meta['prod']['card4l-version'] = '1.0'
    else:
        meta['prod']['card4l-link'] = 'https://ceos.org/ard/files/PFS/NRB/v5.5/CARD4L-PFS_NRB_v5.5.pdf'
        meta['prod']['card4l-version'] = '5.5'
    meta['prod']['compression_type'] = compression
    meta['prod']['compression_zerrors'] = z_err_dict
    meta['prod']['crsEPSG'] = str(prod_meta['epsg'])
    meta['prod']['crsWKT'] = prod_meta['wkt']
    meta['prod']['demAccess'] = DEM_MAP[config['dem_type']]['access']
    meta['prod']['demEGMReference'] = DEM_MAP[config['dem_type']]['egm']
    meta['prod']['demEGMResamplingMethod'] = 'bilinear'
    meta['prod']['demGSD'] = DEM_MAP[config['dem_type']]['gsd']
    meta['prod']['demName'] = config['dem_type'].replace(' II', '')
    meta['prod']['demReference'] = DEM_MAP[config['dem_type']]['ref']
    meta['prod']['demResamplingMethod'] = 'bilinear'
    meta['prod']['demType'] = DEM_MAP[config['dem_type']]['type']
    meta['prod']['doi'] = config['meta']['doi']
    meta['prod']['ellipsoidalHeight'] = None
    meta['prod']['equivalentNumberLooks'] = np.round(calc_enl(tif=ref_tif), 2)
    meta['prod']['geoCorrAccuracyEasternBias'] = None
    meta['prod']['geoCorrAccuracyEasternSTDev'] = None
    meta['prod']['geoCorrAccuracyNorthernBias'] = None
    meta['prod']['geoCorrAccuracyNorthernSTDev'] = None
    meta['prod']['geoCorrAccuracyReference'] = 'https://s1-nrb.readthedocs.io/en/v{}/general/geoaccuracy.html' \
                                               ''.format(S1_NRB.__version__)
    meta['prod']['geoCorrAccuracyType'] = 'slant-range'
    meta['prod']['geoCorrAccuracy_rRMSE'] = geo_corr_accuracy
    meta['prod']['geoCorrAlgorithm'] = 'https://sentinel.esa.int/documents/247904/1653442/' \
                                       'Guide-to-Sentinel-1-Geocoding.pdf'
    meta['prod']['geoCorrResamplingMethod'] = 'bilinear'
    meta['prod']['geom_stac_bbox_native'] = prod_meta['geom']['bbox_native']
    meta['prod']['geom_stac_bbox_4326'] = prod_meta['geom']['bbox']
    meta['prod']['geom_stac_geometry_4326'] = prod_meta['geom']['geometry']
    meta['prod']['geom_xml_center'] = prod_meta['geom']['center']
    meta['prod']['geom_xml_envelope'] = prod_meta['geom']['envelop']
    meta['prod']['griddingConvention'] = 'Military Grid Reference System (MGRS)'
    meta['prod']['griddingConventionURL'] = 'http://www.mgrs-data.org/data/documents/nga_mgrs_doc.pdf'
    meta['prod']['licence'] = config['meta']['licence']
    meta['prod']['mgrsID'] = prod_meta['mgrsID']
    meta['prod']['NRApplied'] = True
    meta['prod']['NRAlgorithm'] = 'https://sentinel.esa.int/documents/247904/2142675/Thermal-Denoising-of-Products-' \
                                  'Generated-by-Sentinel-1-IPF' if meta['prod']['NRApplied'] else None
    meta['prod']['numberOfAcquisitions'] = str(len(src_sid))
    meta['prod']['numBorderPixels'] = prod_meta['nodata_borderpx']
    meta['prod']['numLines'] = str(prod_meta['rows'])
    meta['prod']['numPixelsPerLine'] = str(prod_meta['cols'])
    meta['prod']['pixelCoordinateConvention'] = 'upper-left'
    meta['prod']['processingCenter'] = config['meta']['processing_center']
    meta['prod']['processingMode'] = 'PROTOTYPE'
    meta['prod']['processorName'] = 'S1_NRB'
    meta['prod']['processorVersion'] = S1_NRB.__version__
    meta['prod']['productName'] = 'Ocean Radar Backscatter' if product_type == 'ORB' else 'Normalised Radar Backscatter'
    meta['prod']['productName-short'] = product_type
    meta['prod']['pxSpacingColumn'] = str(prod_meta['res'][0])
    meta['prod']['pxSpacingRow'] = str(prod_meta['res'][1])
    meta['prod']['radiometricAccuracyAbsolute'] = None
    meta['prod']['radiometricAccuracyRelative'] = None
    meta['prod']['radiometricAccuracyReference'] = None
    meta['prod']['rangeNumberOfLooks'] = round(prod_meta['ML_nRgLooks'], 2)
    meta['prod']['RTCAlgorithm'] = 'https://doi.org/10.1109/Tgrs.2011.2120616' \
        if meta['prod']['backscatterMeasurement'] == 'gamma0' or len(ratio_tif) > 0 else None
    meta['prod']['speckleFilterApplied'] = False
    meta['prod']['status'] = 'PLANNED'
    meta['prod']['timeCreated'] = proc_time
    meta['prod']['timeStart'] = start
    meta['prod']['timeStop'] = stop
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
        osv = sid.getOSV(returnMatch=True, osvType=['POE', 'RES'], useLocal=True)
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
        
        inc_vals = dissolve(list(inc.values()))
        pslr, islr = calc_pslr_islr(annotation_dict=src_xml[uid]['annotation'])
        
        def _read_manifest(pattern, attrib=None):
            obj = src_xml[uid]['manifest'].find(pattern, nsmap)
            if attrib is not None:
                return obj.attrib[attrib]
            else:
                return obj.text
        
        # (sorted alphabetically)
        meta['source'][uid] = {}
        meta['source'][uid]['access'] = 'https://scihub.copernicus.eu'
        meta['source'][uid]['acquisitionType'] = 'NOMINAL'
        meta['source'][uid]['ascendingNodeDate'] = _read_manifest('.//s1:ascendingNodeTime')
        meta['source'][uid]['azimuthLookBandwidth'] = az_look_bandwidth
        meta['source'][uid]['azimuthNumberOfLooks'] = az_num_looks
        meta['source'][uid]['azimuthPixelSpacing'] = az_px_spacing
        meta['source'][uid]['azimuthResolution'] = res_az
        meta['source'][uid]['dataGeometry'] = data_geometry
        meta['source'][uid]['datatakeID'] = _read_manifest('.//s1sarl1:missionDataTakeID')
        meta['source'][uid]['doi'] = 'https://sentinel.esa.int/documents/247904/1877131/' \
                                     'Sentinel-1-Product-Specification'
        meta['source'][uid]['faradayMeanRotationAngle'] = None
        meta['source'][uid]['faradayRotationReference'] = None
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
        meta['source'][uid]['orbitDataAccess'] = 'https://scihub.copernicus.eu/gnss'
        meta['source'][uid]['orbitStateVector'] = os.path.basename(osv).replace('.zip', '')
        for osv in list(OSV_MAP.keys()):
            if osv in meta['source'][uid]['orbitStateVector']:
                meta['source'][uid]['orbitDataSource'] = OSV_MAP[osv]
        if len(np_tifs) > 0:
            meta['source'][uid]['perfEstimates'] = calc_performance_estimates(files=np_tifs)
            meta['source'][uid]['perfNoiseEquivalentIntensityType'] = 'sigma0'
        else:
            stats = {stat: None for stat in ['minimum', 'mean', 'maximum']}
            pe = {pol: stats for pol in meta['common']['polarisationChannels']}
            meta['source'][uid]['perfEstimates'] = pe
            meta['source'][uid]['perfNoiseEquivalentIntensityType'] = None
        meta['source'][uid]['perfEquivalentNumberOfLooks'] = enl
        meta['source'][uid]['perfIntegratedSideLobeRatio'] = round(islr, 2)
        meta['source'][uid]['perfPeakSideLobeRatio'] = round(pslr, 2)
        meta['source'][uid]['polCalMatrices'] = None
        fac_org = _read_manifest('.//safe:facility', attrib='organisation')
        fac_name = _read_manifest('.//safe:facility', attrib='name')
        meta['source'][uid]['processingCenter'] = f"{fac_org} {fac_name}".replace(' -', '')
        meta['source'][uid]['processingDate'] = _read_manifest('.//safe:processing', attrib='stop')
        meta['source'][uid]['processingLevel'] = _read_manifest('.//safe:processing', attrib='name')
        meta['source'][uid]['processingMode'] = 'NOMINAL'
        meta['source'][uid]['processorName'] = _read_manifest('.//safe:software', attrib='name')
        meta['source'][uid]['processorVersion'] = _read_manifest('.//safe:software', attrib='version')
        meta['source'][uid]['productType'] = product_type
        meta['source'][uid]['rangeLookBandwidth'] = rg_look_bandwidth
        meta['source'][uid]['rangeNumberOfLooks'] = rg_num_looks
        meta['source'][uid]['rangePixelSpacing'] = rg_px_spacing
        meta['source'][uid]['rangeResolution'] = res_rg
        meta['source'][uid]['sensorCalibration'] = 'https://sentinel.esa.int/web/sentinel/technical-guides/' \
                                                   'sentinel-1-sar/sar-instrument/calibration'
        meta['source'][uid]['status'] = 'ARCHIVED'
        meta['source'][uid]['swaths'] = swaths
        meta['source'][uid]['timeCompletionFromAscendingNode'] = str(float(_read_manifest('.//s1:stopTimeANX')))
        meta['source'][uid]['timeStartFromAscendingNode'] = str(float(_read_manifest('.//s1:startTimeANX')))
        meta['source'][uid]['timeStart'] = datetime.strptime(sid.start, '%Y%m%dT%H%M%S')
        meta['source'][uid]['timeStop'] = datetime.strptime(sid.stop, '%Y%m%dT%H%M%S')
    
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
    
    with _vec_from_srccoords(coord_list=coord_list) as srcvec:
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


def _vec_from_srccoords(coord_list):
    """
    Creates a single :class:`~spatialist.vector.Vector` object from a list of footprint coordinates of source scenes.
    
    Parameters
    ----------
    coord_list: list[list[tuple[float]]]
        List containing for each source scene a list of coordinate pairs as retrieved from the metadata stored in an
        :class:`~pyroSAR.drivers.ID` object.
    
    Returns
    -------
    spatialist.vector.Vector
    """
    if len(coord_list) == 2:
        # determine joined border between footprints
        if math.isclose(coord_list[0][0][0], coord_list[1][3][0], abs_tol=0.1):
            c1 = coord_list[1]
            c2 = coord_list[0]
        elif math.isclose(coord_list[1][0][0], coord_list[0][3][0], abs_tol=0.1):
            c1 = coord_list[0]
            c2 = coord_list[1]
        else:
            raise RuntimeError('Not able to find joined border of source scene footprint coordinates:'
                               '\n{} \n{}'.format(coord_list[0], coord_list[1]))
        
        c1_lat = [c1[0][1], c1[1][1], c1[2][1], c1[3][1]]
        c1_lon = [c1[0][0], c1[1][0], c1[2][0], c1[3][0]]
        c2_lat = [c2[0][1], c2[1][1], c2[2][1], c2[3][1]]
        c2_lon = [c2[0][0], c2[1][0], c2[2][0], c2[3][0]]
        
        wkt = 'POLYGON (({} {},{} {},{} {},{} {},{} {}))'.format(c1_lon[0], c1_lat[0],
                                                                 c1_lon[1], c1_lat[1],
                                                                 c2_lon[2], c2_lat[2],
                                                                 c2_lon[3], c2_lat[3],
                                                                 c1_lon[0], c1_lat[0])
    else:  # len(coord_list) == 1
        c = coord_list[0]
        lat = [c[0][1], c[1][1], c[2][1], c[3][1]]
        lon = [c[0][0], c[1][0], c[2][0], c[3][0]]
        
        wkt = 'POLYGON (({} {},{} {},{} {},{} {},{} {}))'.format(lon[0], lat[0],
                                                                 lon[1], lat[1],
                                                                 lon[2], lat[2],
                                                                 lon[3], lat[3],
                                                                 lon[0], lat[0])
    return wkt2vector(wkt, srs=4326)


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
    files = sid.findfiles(r'^s1[ab].*-[vh]{2}-.*\.xml$')
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
    out['envelop'] = ' '.join(wkt_list)
    
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


def calc_enl(tif, block_size=30, return_arr=False):
    """
    Calculate the Equivalent Number of Looks (ENL) for a linear-scaled backscatter measurement GeoTIFF file. The
    calculation is performed block-wise for the entire image and by default the median ENL value is returned.
    
    Parameters
    ----------
    tif: str
        The path to a linear-scaled backscatter measurement GeoTIFF file.
    block_size: int, optional
        The block size to use for the calculation. Remainder pixels are discarded, if the array dimensions are not
        evenly divisible by the block size. Default is 30, which calculates ENL for 30x30 pixel blocks.
    return_arr: bool, optional
        If True, the calculated ENL array is returned. Default is False.
    
    Returns
    -------
    float or numpy.ndarray
        The median ENL value or array of ENL values if `return_enl_arr` is True.
    
    References
    ----------
    .. [1]  S. N. Anfinsen, A. P. Doulgeris and T. Eltoft,
            "Estimation of the Equivalent Number of Looks in Polarimetric Synthetic Aperture Radar Imagery,"
            in IEEE Transactions on Geoscience and Remote Sensing, vol. 47, no. 11, pp. 3795-3809, Nov. 2009,
            doi: 10.1109/TGRS.2009.2019269.
    """
    with Raster(tif) as ras:
        arr = ras.array()
    arr[np.isinf(arr)] = np.nan
    
    num_blocks_rows = arr.shape[0] // block_size
    num_blocks_cols = arr.shape[1] // block_size
    if num_blocks_rows == 0 or num_blocks_cols == 0:
        raise ValueError("Block size is too large for the input data dimensions.")
    blocks = arr[:num_blocks_rows * block_size,
                 :num_blocks_cols * block_size].reshape(num_blocks_rows, block_size, num_blocks_cols, block_size)
    
    with np.testing.suppress_warnings() as sup:
        sup.filter(RuntimeWarning, "Mean of empty slice")
        _mean = np.nanmean(blocks, axis=(1, 3))
    with np.testing.suppress_warnings() as sup:
        sup.filter(RuntimeWarning, "Degrees of freedom <= 0 for slice")
        _std = np.nanstd(blocks, axis=(1, 3))
    enl = np.divide(_mean ** 2, _std ** 2, out=np.full_like(_mean, fill_value=np.nan), where=_std != 0)
    
    out_arr = np.zeros((num_blocks_rows, num_blocks_cols))
    out_arr[:num_blocks_rows, :num_blocks_cols] = enl
    
    if return_arr:
        return out_arr
    else:
        return np.nanmedian(out_arr)


def calc_geolocation_accuracy(swath_identifier, ei_tif, etad):
    """
    Calculates the radial root mean square error, which is a target requirement of the CARD4L NRB specification
    (Item 4.3). For more information see: https://s1-nrb.readthedocs.io/en/latest/general/geoaccuracy.html.
    Currently only the Copernicus DEM is supported.

    Parameters
    ----------
    swath_identifier: str
        Swath identifier dependent on acquisition mode.
    ei_tif: str
        Path to the annotation GeoTIFF layer 'Ellipsoidal Incident Angle' of the current product.
    etad: bool
        Was the ETAD correction applied?

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
    
    return round(rmse_planar, 2)


def calc_performance_estimates(files):
    """
    Calculates the performance estimates specified in CARD4L NRB 1.6.9 for all noise power images if available.
    
    Parameters
    ----------
    files: list[str]
        List of paths pointing to the noise power images the estimates should be calculated for.
    
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
        out[pol] = {'minimum': round(_min, 2),
                    'maximum': round(_max, 2),
                    'mean': round(_mean, 2)}
    return out


def calc_pslr_islr(annotation_dict):
    """
    Extracts all values for Peak Side Lobe Ratio (PSLR) and Integrated Side Lobe Ratio (ISLR) from the annotation
    metadata of a scene and calculates the mean value for all swaths.
    
    Parameters
    ----------
    annotation_dict: dict
        A dictionary of annotation files in the form: {'swath ID':`lxml.etree._Element` object}
    
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
    pslr = np.nanmean(list(pslr_mean.values()))
    islr = np.nanmean(list(islr_mean.values()))
    
    return pslr, islr


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


def copy_src_meta(target, src_ids):
    """
    Copies the original metadata of the source scenes to the ARD product
    directory.
    
    Parameters
    ----------
    target: str
        A path pointing to the current ARD product directory.
    src_ids: list[pyroSAR.drivers.ID]
        List of :class:`~pyroSAR.drivers.ID` objects of all source scenes that overlap with the current MGRS tile.
    
    Returns
    -------
    None
    """
    for src_id in src_ids:
        source_dir = os.path.join(target, 'source')
        pid = re.match(src_id.pattern, os.path.basename(src_id.file)).group('productIdentifier')
        
        if src_id.scene.endswith('.zip'):
            base = os.path.basename(src_id.file)
            with zipfile.ZipFile(src_id.scene, 'r') as zip_ref:
                zip_ref.extract(member=base + '/manifest.safe', path=source_dir)
                annotation_files = [f for f in zip_ref.namelist() if base + '/annotation' in f]
                zip_ref.extractall(members=annotation_files, path=source_dir)
            os.rename(os.path.join(source_dir, base), os.path.join(source_dir, pid))
        else:
            pid_dir = os.path.join(source_dir, pid)
            os.makedirs(pid_dir, exist_ok=True)
            shutil.copy(src=os.path.join(src_id.scene, 'manifest.safe'),
                        dst=os.path.join(pid_dir, 'manifest.safe'))
            shutil.copytree(src=os.path.join(src_id.scene, 'annotation'),
                            dst=os.path.join(pid_dir, 'annotation'))
