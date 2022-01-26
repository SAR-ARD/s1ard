import os
import re
from lxml import etree
from datetime import datetime
from pyroSAR import identify
from spatialist import Raster
from spatialist.ancillary import finder
from spatialist.vector import wkt2vector
from spatialist.raster import rasterize
from math import isclose
import numpy as np

from S1_NRB.metadata.mapping import NRB_PATTERN, RES_MAP, ORB_MAP


def get_prod_meta(product_id, tif, sources):
    """
    Returns a metadata dictionary, which is generated from the ID of a product scene using a regular expression pattern
    and from a measurement GeoTIFF file of the same product scene using spatialist's Raster class.
    
    Parameters
    ----------
    product_id: str
        The product ID (filename) of the product scene.
    tif: str
        The paths to a measurement GeoTIFF file of the product scene.
    sources: list[str]
        A list of paths pointing to the source scenes of the product.
    
    Returns
    -------
    dict
        A dictionary containing metadata for the product scene.
    """
    
    out = re.match(re.compile(NRB_PATTERN), product_id).groupdict()
    coord_list = [identify(src).meta['coordinates'] for src in sources]
    
    with vec_from_srccoords(coord_list=coord_list) as srcvec:
        with Raster(tif) as ras:
            vec = ras.bbox()
            srs = vec.srs
            out['extent'] = vec.extent
            out['wkt'] = srs.ExportToWkt()
            out['epsg'] = vec.getProjection(type='epsg')
            out['rows'] = ras.rows
            out['cols'] = ras.cols
            out['res'] = ras.res
            geo = ras.geo
            out['transform'] = [geo['xres'], geo['rotation_x'], geo['xmin'],
                                geo['rotation_y'], geo['yres'], geo['ymax']]
            
            vec.reproject(4326)
            feat = vec.getFeatureByIndex(0)
            geom = feat.GetGeometryRef()
            point = geom.Centroid()
            out['wkt_pt'] = point.ExportToWkt()
            out['wkt_env'] = vec.convert2wkt(set3D=False)[0]
            out['extent_4326'] = vec.extent
            
            # Calculate number of nodata border pixels based on source scene(s) footprint
            ras_srcvec = rasterize(vectorobject=srcvec, reference=ras, burn_values=[1])
            arr_srcvec = ras_srcvec.array()
            out['nodata_borderpx'] = np.count_nonzero(np.isnan(arr_srcvec))
    
    return out


def get_uid_sid(filepath):
    """
    Returns the unique identifier of a Sentinel-1 scene and a pyroSAR metadata handler generated using 
    `pyroSAR.drivers.identify`
    
    Parameters
    ----------
    filepath: str
        Filepath pointing to a Sentinel-1 scene.
    
    Returns
    -------
    uid: str
        The last four characters of a Sentinel-1 filename, which is the unique identifier of the scene.
    sid: pyroSAR.drivers.ID subclass object
        A pyroSAR metadata handler for the scene generated with `pyroSAR.drivers.identify`.
    """
    
    uid = os.path.basename(filepath).split('.')[0][-4:]
    sid = identify(filepath)
    
    return uid, sid


def etree_from_sid(sid):
    """
    Uses a pyroSAR metadata handler to get the parsed manifest and annotation XML data as a dictionary.
    
    Parameters
    ----------
    sid:  pyroSAR.drivers.ID subclass object
        A pyroSAR metadata handler generated with `pyroSAR.drivers.identify`.
    
    Returns
    -------
    dict
        A dictionary containing the parsed etree.ElementTree objects for the manifest and annotation XML files.
    """
    
    with sid.getFileObj(sid.findfiles('manifest.safe')[0]) as input_man:
        annotation_files = sid.findfiles(r'^s1[ab].*-vh-.*\.xml$')
        swaths = [re.search('-iw[1-3]|-ew[1-5]|-s[1-6]]', a).group().replace('-', '') for a in annotation_files]
        annotation_dict = {}
        for s, a in zip(swaths, annotation_files):
            annotation_dict[s.upper()] = etree.fromstring(sid.getFileObj(a).getvalue())
        
        return {'manifest': etree.fromstring(input_man.getvalue()),
                'annotation': annotation_dict}


def convert_id_coordinates(coords, stac=False):
    """
    Converts a list of coordinate pairs as retrieved by pyroSAR's identify function to either envelop and center for
    usage in the XML metadata files or bbox and geometry for usage in STAC metadata files. The latter is returned if the
    optional parameter 'stac' is set to True, else the former is returned.
    
    Parameters
    ----------
    coords: list[tuple(float, float)]
        List of coordinate pairs as retrieved by `pyroSAR.drivers.identify` from source scenes
    stac: bool, optional
        If set to True, bbox and geometry are returned for usage in STAC Items. If set to False (default) envelop and
        center are returned for usage in XML metadata files.
    
    Returns
    -------
    envelop: str
        Acquisition footprint coordinates for the XML field 'multiExtentOf' in 'eop:Footprint'
    center: str
        Acquisition center coordinates for the XML field 'centerOf' in 'eop:Footprint'
    
    Notes
    -------
    If stac=True the following parameters are returned instead of envelop and center:
    
    bbox: list[float]
        Acquisition bounding box for usage in STAC Items. Formatted in accordance with RFC 7946, section 5:
        https://datatracker.ietf.org/doc/html/rfc7946#section-5
    geometry: GeoJSON Geometry Object
        Acquisition footprint geometry for usage in STAC Items. Formatted in accordance with RFC 7946, section 3.1.:
        https://datatracker.ietf.org/doc/html/rfc7946#section-3.1
    """
    c = coords
    
    lat = [c[0][1], c[1][1], c[2][1], c[3][1]]
    lon = [c[0][0], c[1][0], c[2][0], c[3][0]]
    envelop = '{} {},{} {},{} {},{} {},{} {}'.format(lon[0], lat[0],
                                                     lon[1], lat[1],
                                                     lon[2], lat[2],
                                                     lon[3], lat[3],
                                                     lon[0], lat[0])
    
    lon_c = (max(lon) + min(lon)) / 2
    lat_c = (max(lat) + min(lat)) / 2
    center = '{} {}'.format(lon_c, lat_c)
    
    if stac:
        bbox = [min(lon), min(lat), max(lon), max(lat)]
        geometry = {'type': 'Polygon', 'coordinates': (((lon[0], lat[0]),
                                                        (lon[1], lat[1]),
                                                        (lon[2], lat[2]),
                                                        (lon[3], lat[3]),
                                                        (lon[0], lat[0])),)}
        return bbox, geometry
    else:
        return envelop, center


def convert_spatialist_extent(extent):
    """
    Converts the extent of a spatialist vector object to bbox and geometry as required for the usage in STAC Items:
    https://github.com/radiantearth/stac-spec/blob/master/item-spec/item-spec.md#item-fields
    
    Parameters
    ----------
    extent: dict
        The extent of a vector object as returned by `spatialist.vector.Vector.extent`
    
    Returns
    -------
    bbox: list[float]
        Formatted in accordance with RFC 7946, section 5: https://datatracker.ietf.org/doc/html/rfc7946#section-5
    geometry: GeoJSON Geometry Object
        Formatted in accordance with RFC 7946, section 3.1.: https://datatracker.ietf.org/doc/html/rfc7946#section-3.1
    """
    xmin = extent['xmin']
    xmax = extent['xmax']
    ymin = extent['ymin']
    ymax = extent['ymax']
    
    bbox = [xmin, ymin, xmax, ymax]
    geometry = {'type': 'Polygon', 'coordinates': (((xmin, ymin),
                                                    (xmin, ymax),
                                                    (xmax, ymax),
                                                    (xmax, ymin),
                                                    (xmin, ymin)),)}
    
    return bbox, geometry


def vec_from_srccoords(coord_list):
    """
    Creates a single `spatialist.vector.Vector` object from a list of footprint coordinates of source scenes.
    
    Parameters
    ----------
    coord_list: list[list[tuple(float, float)]]
        List containing (for n source scenes) a list of coordinate pairs as retrieved by `pyroSAR.drivers.identify`
    
    Returns
    -------
    `spatialist.vector.Vector` object
    """
    if len(coord_list) == 2:
        if isclose(coord_list[0][0][0], coord_list[1][3][0], abs_tol=0.1):
            c1 = coord_list[1]
            c2 = coord_list[0]
        elif isclose(coord_list[1][0][0], coord_list[0][3][0], abs_tol=0.1):
            c1 = coord_list[0]
            c2 = coord_list[1]
        else:
            RuntimeError('not able to find joined border of source scene footprints '
                         '\n{} and \n{}'.format(coord_list[0], coord_list[1]))
        
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


def find_in_annotation(annotation_dict, pattern, single=False, out_type=None):
    """
    Find a pattern in each annotation file provided and returns a list of results converted in
    
    Parameters
    ----------
    annotation_dict: dict
        A dict of annotation files in the form: {'swath ID': lxml.etree._Element object}
    pattern: str
        The pattern to search for in each annotation file.
    single: bool
        If True, the results found in each annotation file are expected to be the same and therefore only a single
        value will be returned instead of a dict. Default is False.
    out_type: str, optional
        Output type to convert the results to. Can be one of the following:
        str (default)
        float
        int
    
    Returns
    -------
    out: dict
        A dict of the results containing a list for each of the annotation files.
        I.e., {'swath ID': list[str, float or int]}
    """
    if out_type is None:
        out_type = 'str'
    
    out = {}
    for s, a in annotation_dict.items():
        out[s] = [x.text for x in a.findall(pattern)]
        if len(out[s]) == 1:
            out[s] = out[s][0]
    
    if out_type != 'str':
        for k in list(out.keys()):
            if isinstance(out[k], list):
                if out_type == 'float':
                    out[k] = [float(x) for x in out[k]]
                elif out_type == 'int':
                    out[k] = [int(x) for x in out[k]]
            else:
                if out_type == 'float':
                    out[k] = float(out[k])
                elif out_type == 'int':
                    out[k] = int(out[k])
    
    if single:
        test_val = list(out.values())[0]
        for k in out:
            if out[k] != test_val:
                raise RuntimeError('Search result for pattern "{}" expected to be the same in all annotation '
                                   'files.'.format(pattern))
        if out_type == 'float':
            return float(test_val)
        elif out_type == 'int':
            return int(test_val)
        else:
            return str(test_val)
    else:
        return out


def meta_dict(target, src_scenes, src_files, dem_name, proc_time):
    """
    Creates a dictionary containing metadata for a product scene, as well as its source scenes. The dictionary can then
    be utilized by `metadata.xmlparser` and `metadata.stacparser` to generate XML and STAC JSON metadata files,
    respectively.
    
    Parameters
    ----------
    target: str
        A path pointing to the root directory of a product scene.
    src_scenes: list[str]
        A list of paths pointing to the source scenes of the product.
    src_files: list[str]
        A list of paths pointing to the SNAP processed datasets of the product.
    dem_name: str
        Name of the DEM used for processing.
    proc_time: datetime.datetime
        The datetime object used to generate the unique product identifier from.
    
    Returns
    -------
    meta: dict
        A dictionary containing an extensive collection of metadata for product as well as source scenes.
    """
    
    meta = {'prod': {},
            'source': {},
            'common': {}}
    
    product_id = os.path.basename(target)
    tif = finder(target, ['vh-g-lin.tif$'], regex=True)[0]
    prod_meta = get_prod_meta(product_id=product_id, tif=tif, sources=src_scenes)
    
    src_sid = {}
    src_xml = {}
    for i in range(len(src_scenes)):
        uid, sid = get_uid_sid(filepath=src_scenes[i])
        src_sid[uid] = sid
        src_xml[uid] = etree_from_sid(sid=sid)
    
    src0 = list(src_sid.keys())[0]  # first key/first file
    sid0 = src_sid[src0]
    manifest0 = src_xml[src0]['manifest']
    nsmap0 = src_xml[src0]['manifest'].nsmap
    
    stac_bbox_4326, stac_geometry_4326 = convert_spatialist_extent(extent=prod_meta['extent_4326'])
    stac_bbox_native = convert_spatialist_extent(extent=prod_meta['extent'])[0]
    
    dem_map = {'GETASSE30': {'url': 'https://step.esa.int/auxdata/dem/GETASSE30',
                             'ref': 'https://seadas.gsfc.nasa.gov/help-8.1.0/desktop/GETASSE30ElevationModel.html',
                             'type': 'elevation',
                             'egm': 'https://apps.dtic.mil/sti/citations/ADA166519'},
               'Copernicus 10m EEA DEM': {'url': 'ftps://cdsdata.copernicus.eu/DEM-datasets/COP-DEM_EEA-10-DGED/',
                                          'ref': 'https://spacedata.copernicus.eu/web/cscda/dataset-details?articleId=394198',
                                          'type': 'surface',
                                          'egm': 'https://doi.org/10.1029/2011JB008916'},
               'Copernicus 30m Global DEM': {'url': 'https://registry.opendata.aws/copernicus-dem/',
                                             'ref': 'https://spacedata.copernicus.eu/web/cscda/dataset-details?articleId=394198',
                                             'type': 'surface',
                                             'egm': 'https://doi.org/10.1029/2011JB008916'}}
    dem_url = dem_map[dem_name]['url']
    dem_ref = dem_map[dem_name]['ref']
    dem_type = dem_map[dem_name]['type']
    dem_egm = dem_map[dem_name]['egm']
    
    # Common metadata (sorted alphabetically)
    meta['common']['antennaLookDirection'] = 'RIGHT'
    meta['common']['constellation'] = 'sentinel-1'
    meta['common']['instrumentShortName'] = 'C-SAR'
    meta['common']['operationalMode'] = prod_meta['mode']
    meta['common']['orbitMeanAltitude'] = '{:.2e}'.format(693000)
    meta['common']['orbit'] = {'A': 'ascending', 'D': 'descending'}[sid0.orbit]
    meta['common']['orbitNumbers_abs'] = sid0.meta['orbitNumbers_abs']
    meta['common']['orbitNumbers_rel'] = sid0.meta['orbitNumbers_rel']
    meta['common']['orbitNumber_start'] = str(meta['common']['orbitNumbers_abs']['start'])
    meta['common']['orbitNumber_stop'] = str(meta['common']['orbitNumbers_abs']['stop'])
    meta['common']['platformIdentifier'] = {'S1A': '1A', 'S1B': '1B'}[sid0.sensor]
    meta['common']['platformShortName'] = 'Sentinel'
    meta['common']['platformFullname'] = '{}-{}'.format(meta['common']['platformShortName'].lower(),
                                                        meta['common']['platformIdentifier'].lower())
    meta['common']['platformReference'] = {'sentinel-1a': 'http://database.eohandbook.com/database/missionsummary.aspx?missionID=575',
                                           'sentinel-1b': 'http://database.eohandbook.com/database/missionsummary.aspx?missionID=576'}[meta['common']['platformFullname']]
    meta['common']['polarisationChannels'] = sid0.polarizations
    meta['common']['polarisationMode'] = prod_meta['pols']
    meta['common']['radarBand'] = 'C'
    meta['common']['radarCenterFreq'] = '{:.3e}'.format(5405000000)
    meta['common']['sensorType'] = 'RADAR'
    
    # Product metadata (sorted alphabetically)
    meta['prod']['access'] = None
    meta['prod']['ancillaryData1'] = None
    meta['prod']['acquisitionType'] = 'NOMINAL'
    meta['prod']['ascendingNodeDate'] = manifest0.find('.//s1:ascendingNodeTime', nsmap0).text
    meta['prod']['backscatterConvention'] = 'linear power'
    meta['prod']['backscatterConversionEq'] = '10*log10(DN)'
    meta['prod']['backscatterMeasurement'] = 'gamma0'
    meta['prod']['card4l-link'] = 'https://ceos.org/ard/files/PFS/NRB/v5.5/CARD4L-PFS_NRB_v5.5.pdf'
    meta['prod']['card4l-name'] = 'NRB'
    meta['prod']['card4l-version'] = '5.5'
    meta['prod']['crsEPSG'] = str(prod_meta['epsg'])
    meta['prod']['crsWKT'] = prod_meta['wkt']
    meta['prod']['datatakeID'] = manifest0.find('.//s1sarl1:missionDataTakeID', nsmap0).text
    meta['prod']['demEgmReference'] = dem_egm
    meta['prod']['demEgmResamplingMethod'] = 'bilinear'
    meta['prod']['demName'] = dem_name
    meta['prod']['demReference'] = dem_ref
    meta['prod']['demResamplingMethod'] = 'bilinear'
    meta['prod']['demType'] = dem_type
    meta['prod']['demURL'] = dem_url
    meta['prod']['doi'] = None
    meta['prod']['ellipsoidalHeight'] = None
    meta['prod']['fileBitsPerSample'] = '32'
    meta['prod']['fileByteOrder'] = 'little-endian'
    meta['prod']['fileDataType'] = 'float'
    meta['prod']['fileFormat'] = 'COG'
    meta['prod']['filterApplied'] = False
    meta['prod']['filterType'] = None
    meta['prod']['filterWindowSizeCol'] = None
    meta['prod']['filterWindowSizeLine'] = None
    meta['prod']['geoCorrAccuracyEasternBias'] = None
    meta['prod']['geoCorrAccuracyEasternSTDev'] = None
    meta['prod']['geoCorrAccuracyNorthernBias'] = None
    meta['prod']['geoCorrAccuracyNorthernSTDev'] = None
    meta['prod']['geoCorrAccuracy_rRMSE'] = None
    meta['prod']['geoCorrAccuracyReference'] = 'https://www.mdpi.com/2072-4292/9/6/607'
    meta['prod']['geoCorrAccuracyType'] = 'slant-range'
    meta['prod']['geoCorrAlgorithm'] = 'https://sentinel.esa.int/documents/247904/1653442/Guide-to-Sentinel-1-Geocoding.pdf'
    meta['prod']['geoCorrResamplingMethod'] = 'bilinear'
    meta['prod']['geom_stac_bbox_native'] = stac_bbox_native
    meta['prod']['geom_stac_bbox_4326'] = stac_bbox_4326
    meta['prod']['geom_stac_geometry_4326'] = stac_geometry_4326
    meta['prod']['geom_xml_center'] = re.search(r'\(([-*0-9 .,]+)\)', prod_meta['wkt_pt']).group(1)
    meta['prod']['geom_xml_envelope'] = re.search(r'\(([-*0-9 .,]+)\)', prod_meta['wkt_env']).group(1)
    meta['prod']['griddingConventionURL'] = 'http://www.mgrs-data.org/data/documents/nga_mgrs_doc.pdf'
    meta['prod']['griddingConvention'] = 'Military Grid Reference System (MGRS)'
    meta['prod']['licence'] = None
    meta['prod']['majorCycleID'] = str(sid0.meta['cycleNumber'])
    meta['prod']['noiseRemovalApplied'] = True
    meta['prod']['noiseRemovalAlgorithm'] = 'https://doi.org/10.1109/tgrs.2018.2889381'
    meta['prod']['numberLines'] = str(prod_meta['rows'])
    meta['prod']['numberOfAcquisitions'] = str(len(src_scenes))
    meta['prod']['numBorderPixels'] = prod_meta['nodata_borderpx']
    meta['prod']['numPixelsPerLine'] = str(prod_meta['cols'])
    meta['prod']['pixelCoordinateConvention'] = 'pixel ULC'
    meta['prod']['processingCenter'] = 'FSU'
    meta['prod']['processingLevel'] = 'L2'
    meta['prod']['processingMode'] = 'PROTOTYPE'
    meta['prod']['processorName'] = 'S1_NRB'
    meta['prod']['processorVersion'] = '0.1'
    meta['prod']['productName'] = 'NORMALISED RADAR BACKSCATTER'
    meta['prod']['pxSpacingColumn'] = str(prod_meta['res'][0])
    meta['prod']['pxSpacingRow'] = str(prod_meta['res'][1])
    meta['prod']['radiometricAccuracyAbsolute'] = None
    meta['prod']['radiometricAccuracyRelative'] = None
    meta['prod']['radiometricAccuracyReference'] = None
    meta['prod']['RTCAlgorithm'] = 'https://doi.org/10.1109/Tgrs.2011.2120616'
    meta['prod']['status'] = 'PROTOTYPE'
    meta['prod']['timeCreated'] = proc_time
    meta['prod']['timeCompletionFromAscendingNode'] = str(float(manifest0.find('.//s1:stopTimeANX', nsmap0).text))
    meta['prod']['timeStartFromAscendingNode'] = str(float(manifest0.find('.//s1:startTimeANX', nsmap0).text))
    meta['prod']['timeStart'] = datetime.strptime(prod_meta['start'], '%Y%m%dT%H%M%S')
    meta['prod']['timeStop'] = datetime.strptime(prod_meta['stop'], '%Y%m%dT%H%M%S')
    meta['prod']['transform'] = prod_meta['transform']
    meta['prod']['wrsLongitudeGrid'] = str(meta['common']['orbitNumbers_rel']['start'])
    
    # Source metadata
    for uid in list(src_sid.keys()):
        nsmap = src_xml[uid]['manifest'].nsmap
        
        osv = src_sid[uid].getOSV(returnMatch=True, osvType=['POE', 'RES'], useLocal=True)
        coords = src_sid[uid].meta['coordinates']
        xml_envelop, xml_center = convert_id_coordinates(coords=coords)
        stac_bbox, stac_geometry = convert_id_coordinates(coords=coords, stac=True)
        swaths = list(src_xml[uid]['annotation'].keys())
        
        # (sorted alphabetically)
        meta['source'][uid] = {}
        meta['source'][uid]['access'] = 'https://scihub.copernicus.eu'
        meta['source'][uid]['acquisitionType'] = 'NOMINAL'
        meta['source'][uid]['azimuthLookBandwidth'] = find_in_annotation(annotation_dict=src_xml[uid]['annotation'],
                                                                         pattern='.//azimuthProcessing/lookBandwidth',
                                                                         out_type='float')
        meta['source'][uid]['azimuthNumberOfLooks'] = find_in_annotation(annotation_dict=src_xml[uid]['annotation'],
                                                                         pattern='.//azimuthProcessing/numberOfLooks',
                                                                         single=True)
        try:
            meta['source'][uid]['azimuthPixelSpacing'] = find_in_annotation(annotation_dict=src_xml[uid]['annotation'],
                                                                            pattern='.//azimuthPixelSpacing',
                                                                            single=True)
        except RuntimeError:
            tmp_out = find_in_annotation(annotation_dict=src_xml[uid]['annotation'],
                                         pattern='.//azimuthPixelSpacing',
                                         single=False, out_type='float')
            meta['source'][uid]['azimuthPixelSpacing'] = str(sum(list(tmp_out.values())) / len(list(tmp_out.values())))
        meta['source'][uid]['azimuthResolution'] = RES_MAP[meta['common']['operationalMode']]['azimuthResolution']
        meta['source'][uid]['dataGeometry'] = 'slant range'
        meta['source'][uid]['doi'] = 'https://sentinel.esa.int/documents/247904/1877131/Sentinel-1-Product-Specification'
        meta['source'][uid]['faradayMeanRotationAngle'] = None
        meta['source'][uid]['faradayRotationReference'] = None
        meta['source'][uid]['filename'] = src_sid[uid].file
        meta['source'][uid]['geom_stac_bbox_4326'] = stac_bbox
        meta['source'][uid]['geom_stac_geometry_4326'] = stac_geometry
        meta['source'][uid]['geom_xml_center'] = xml_center
        meta['source'][uid]['geom_xml_envelop'] = xml_envelop
        inc = find_in_annotation(annotation_dict=src_xml[uid]['annotation'],
                                 pattern='.//geolocationGridPoint/incidenceAngle', out_type='float')
        meta['source'][uid]['incidenceAngleMax'] = np.max(list(inc.values()))
        meta['source'][uid]['incidenceAngleMin'] = np.min(list(inc.values()))
        meta['source'][uid]['incidenceAngleMidSwath'] = np.max(list(inc.values())) - \
                                                        ((np.max(list(inc.values())) - np.min(list(inc.values()))) / 2)
        meta['source'][uid]['instrumentAzimuthAngle'] = str(src_sid[uid].meta['heading'])
        meta['source'][uid]['ionosphereIndicator'] = None
        meta['source'][uid]['lutApplied'] = find_in_annotation(annotation_dict=src_xml[uid]['annotation'],
                                                               pattern='.//applicationLutId',
                                                               single=True)
        meta['source'][uid]['orbitStateVector'] = os.path.basename(osv).replace('.zip', '')
        for orb in list(ORB_MAP.keys()):
            if orb in meta['source'][uid]['orbitStateVector']:
                meta['source'][uid]['orbitDataSource'] = ORB_MAP[orb]
        meta['source'][uid]['orbitDataAccess'] = 'https://scihub.copernicus.eu/gnss'
        meta['source'][uid]['perfEstimatesMax'] = None
        meta['source'][uid]['perfEstimatesMean'] = None
        meta['source'][uid]['perfEstimatesMin'] = None
        meta['source'][uid]['perfEquivalentNumberOfLooks'] = None
        #meta['source'][uid]['perfIndicatorsURL'] = 'https://sentinel.esa.int/documents/247904/2142675/Thermal-Denoising-of-Products-Generated-by-Sentinel-1-IPF'
        meta['source'][uid]['perfIntegratedSideLobeRatio'] = None
        meta['source'][uid]['perfNoiseEquivalentIntensityType'] = 'gamma0'
        meta['source'][uid]['perfPeakSideLobeRatio'] = None
        meta['source'][uid]['polCalMatrices'] = None
        meta['source'][uid]['processingCenter'] = f"{src_xml[uid]['manifest'].find('.//safe:facility', nsmap).attrib['organisation']}, " \
                                                  f"{src_xml[uid]['manifest'].find('.//safe:facility', nsmap).attrib['name']}, " \
                                                  f"{src_xml[uid]['manifest'].find('.//safe:facility', nsmap).attrib['site']}, " \
                                                  f"{src_xml[uid]['manifest'].find('.//safe:facility', nsmap).attrib['country']}"
        meta['source'][uid]['processingDate'] = src_xml[uid]['manifest'].find('.//safe:processing', nsmap).attrib['stop']
        meta['source'][uid]['processingLevel'] = src_xml[uid]['manifest'].find('.//safe:processing', nsmap).attrib['name']
        meta['source'][uid]['processorName'] = src_xml[uid]['manifest'].find('.//safe:software', nsmap).attrib['name']
        meta['source'][uid]['processorVersion'] = src_xml[uid]['manifest'].find('.//safe:software', nsmap).attrib['version']
        meta['source'][uid]['productType'] = src_sid[uid].meta['product']
        meta['source'][uid]['rangeLookBandwidth'] = find_in_annotation(annotation_dict=src_xml[uid]['annotation'],
                                                                       pattern='.//rangeProcessing/lookBandwidth',
                                                                       out_type='float')
        meta['source'][uid]['rangeNumberOfLooks'] = find_in_annotation(annotation_dict=src_xml[uid]['annotation'],
                                                                       pattern='.//rangeProcessing/numberOfLooks',
                                                                       single=True)
        try:
            meta['source'][uid]['rangePixelSpacing'] = find_in_annotation(annotation_dict=src_xml[uid]['annotation'],
                                                                          pattern='.//rangePixelSpacing',
                                                                          single=True)
        except RuntimeError:
            tmp_out = find_in_annotation(annotation_dict=src_xml[uid]['annotation'],
                                         pattern='.//rangePixelSpacing',
                                         single=False, out_type='float')
            meta['source'][uid]['rangePixelSpacing'] = str(sum(list(tmp_out.values())) / len(list(tmp_out.values())))
        meta['source'][uid]['rangeResolution'] = RES_MAP[meta['common']['operationalMode']]['rangeResolution']
        meta['source'][uid]['sensorCalibration'] = 'https://sentinel.esa.int/web/sentinel/technical-guides/sentinel-1-sar/sar-instrument/calibration'
        meta['source'][uid]['status'] = 'ARCHIVED'
        meta['source'][uid]['timeStart'] = datetime.strptime(src_sid[uid].start, '%Y%m%dT%H%M%S')
        meta['source'][uid]['timeStop'] = datetime.strptime(src_sid[uid].stop, '%Y%m%dT%H%M%S')
        meta['source'][uid]['swathIdentifier'] = re.search('_IW|EW|S[1-6]{1}_', os.path.basename(src_sid[uid].file)).group().replace('_', '')
        meta['source'][uid]['swaths'] = swaths
    
    # return meta, m_sid, m_src
    return meta
