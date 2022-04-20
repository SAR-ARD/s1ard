import os
import re
from copy import deepcopy
from lxml import etree
from datetime import datetime
from spatialist import Raster
from statistics import mean
from S1_NRB.metadata.mapping import SAMPLE_MAP, NS_MAP
from S1_NRB.metadata.extract import get_header_size


def _nsc(text, nsmap):
    ns, key = text.split(':')
    return '{{{0}}}{1}'.format(nsmap[ns], key)


def om_time(root, nsmap, scene_id, time_start, time_stop):
    """
    Creates the `om:phenomenonTime` and `om:resultTime` XML elements.
    
    Parameters
    ----------
    root: etree.Element
        Root XML element.
    nsmap: dict
        Dictionary listing abbreviation (key) and URI (value) of all necessary XML namespaces.
    scene_id: str
        Scene basename.
    time_start: str
        Start time of the scene acquisition.
    time_stop: str
        Stop time of the acquisition.
    
    Returns
    -------
    None
    """
    phenomenonTime = etree.SubElement(root, _nsc('om:phenomenonTime', nsmap))
    timePeriod = etree.SubElement(phenomenonTime, _nsc('gml:TimePeriod', nsmap),
                                  attrib={_nsc('gml:id', nsmap): scene_id + '_2'})
    beginPosition = etree.SubElement(timePeriod, _nsc('gml:beginPosition', nsmap))
    beginPosition.text = time_start
    endPosition = etree.SubElement(timePeriod, _nsc('gml:endPosition', nsmap))
    endPosition.text = time_stop
    
    resultTime = etree.SubElement(root, _nsc('om:resultTime', nsmap))
    timeInstant = etree.SubElement(resultTime, _nsc('gml:TimeInstant', nsmap),
                                   attrib={_nsc('gml:id', nsmap): scene_id + '_3'})
    timePosition = etree.SubElement(timeInstant, _nsc('gml:timePosition', nsmap))
    timePosition.text = time_stop


def om_procedure(root, nsmap, scene_id, meta, uid=None, prod=True):
    """
    Creates the `om:procedure/eop:EarthObservationEquipment` XML elements and all relevant subelements for source and
    product metadata. Differences between source and product are controlled using the `prod=[True|False]` switch.
    
    Parameters
    ----------
    root: etree.Element
        Root XML element.
    nsmap: dict
        Dictionary listing abbreviation (key) and URI (value) of all necessary XML namespaces.
    scene_id: str
        Scene basename.
    meta: dict
        Metadata dictionary generated with `metadata.extract.meta_dict`
    uid: str, optional
        Unique identifier of a source SLC scene.
    prod: bool, optional
        Return XML subelements for further usage in `product_xml` parsing function? Default is True. If False, the
        XML subelements for further usage in the `source_xml` parsing function will be returned.
    
    Returns
    -------
    None
    """
    procedure = etree.SubElement(root, _nsc('om:procedure', nsmap))
    earthObservationEquipment = etree.SubElement(procedure, _nsc('eop:EarthObservationEquipment', nsmap),
                                                 attrib={_nsc('gml:id', nsmap): scene_id + '_4'})
    
    # eop:platform
    platform0 = etree.SubElement(earthObservationEquipment, _nsc('eop:platform', nsmap))
    if prod:
        platform1 = etree.SubElement(platform0, _nsc('eop:Platform', nsmap))
    else:
        platform1 = etree.SubElement(platform0, _nsc('nrb:Platform', nsmap))
    shortName = etree.SubElement(platform1, _nsc('eop:shortName', nsmap))
    shortName.text = meta['common']['platformShortName'].upper()
    serialIdentifier = etree.SubElement(platform1, _nsc('eop:serialIdentifier', nsmap))
    serialIdentifier.text = meta['common']['platformIdentifier']
    if not prod:
        satReference = etree.SubElement(platform1, _nsc('nrb:satelliteReference', nsmap),
                                        attrib={_nsc('xlink:href', nsmap): meta['common']['platformReference']})
    
    # eop:instrument
    instrument0 = etree.SubElement(earthObservationEquipment, _nsc('eop:instrument', nsmap))
    instrument1 = etree.SubElement(instrument0, _nsc('eop:Instrument', nsmap))
    shortName = etree.SubElement(instrument1, _nsc('eop:shortName', nsmap))
    shortName.text = meta['common']['instrumentShortName']
    
    # eop:sensor
    sensor0 = etree.SubElement(earthObservationEquipment, _nsc('eop:sensor', nsmap))
    sensor1 = etree.SubElement(sensor0, _nsc('nrb:Sensor', nsmap))
    sensorType = etree.SubElement(sensor1, _nsc('eop:sensorType', nsmap))
    sensorType.text = meta['common']['sensorType']
    operationalMode = etree.SubElement(sensor1, _nsc('eop:operationalMode', nsmap),
                                       attrib={'codeSpace': 'urn:esa:eop:C-SAR:operationalMode'})
    operationalMode.text = meta['common']['operationalMode']
    swathIdentifier = etree.SubElement(sensor1, _nsc('eop:swathIdentifier', nsmap),
                                       attrib={'codeSpace': 'urn:esa:eop:C-SAR:swathIdentifier'})
    swathIdentifier.text = meta['common']['swathIdentifier']
    radarBand = etree.SubElement(sensor1, _nsc('nrb:radarBand', nsmap))
    radarBand.text = meta['common']['radarBand']
    if not prod:
        radarCenterFreq = etree.SubElement(sensor1, _nsc('nrb:radarCenterFrequency', nsmap),
                                           attrib={'uom': 'Hz'})
        radarCenterFreq.text = '{:.3e}'.format(meta['common']['radarCenterFreq'])
        sensorCalibration = etree.SubElement(sensor1, _nsc('nrb:sensorCalibration', nsmap),
                                             attrib={
                                                 _nsc('xlink:href', nsmap): meta['source'][uid]['sensorCalibration']})
    
    # eop:acquisitionParameters
    acquisitionParameters = etree.SubElement(earthObservationEquipment, _nsc('eop:acquisitionParameters', nsmap))
    acquisition = etree.SubElement(acquisitionParameters, _nsc('nrb:Acquisition', nsmap))
    orbitNumber = etree.SubElement(acquisition, _nsc('eop:orbitNumber', nsmap))
    orbitNumber.text = meta['common']['orbitNumber']
    orbitDirection = etree.SubElement(acquisition, _nsc('eop:orbitDirection', nsmap))
    orbitDirection.text = meta['common']['orbitDirection'].upper()
    wrsLongitudeGrid = etree.SubElement(acquisition, _nsc('eop:wrsLongitudeGrid', nsmap),
                                        attrib={'codeSpace': 'urn:esa:eop:Sentinel1:relativeOrbits'})
    wrsLongitudeGrid.text = meta['common']['wrsLongitudeGrid']
    if not prod:
        ascendingNodeDate = etree.SubElement(acquisition, _nsc('eop:ascendingNodeDate', nsmap))
        ascendingNodeDate.text = meta['source'][uid]['ascendingNodeDate']
        startTimeFromAscendingNode = etree.SubElement(acquisition, _nsc('eop:startTimeFromAscendingNode', nsmap),
                                                      attrib={'uom': 'ms'})
        startTimeFromAscendingNode.text = meta['source'][uid]['timeStartFromAscendingNode']
        completionTimeFromAscendingNode = etree.SubElement(acquisition,
                                                           _nsc('eop:completionTimeFromAscendingNode', nsmap),
                                                           attrib={'uom': 'ms'})
        completionTimeFromAscendingNode.text = meta['source'][uid]['timeCompletionFromAscendingNode']
        instrumentAzimuthAngle = etree.SubElement(acquisition, _nsc('eop:instrumentAzimuthAngle', nsmap),
                                                  attrib={'uom': 'deg'})
        instrumentAzimuthAngle.text = meta['source'][uid]['instrumentAzimuthAngle']
    polarisationMode = etree.SubElement(acquisition, _nsc('sar:polarisationMode', nsmap))
    polarisationMode.text = meta['common']['polarisationMode']
    polarisationChannels = etree.SubElement(acquisition, _nsc('sar:polarisationChannels', nsmap))
    polarisationChannels.text = ', '.join(meta['common']['polarisationChannels'])
    if prod:
        numberOfAcquisitions = etree.SubElement(acquisition, _nsc('nrb:numberOfAcquisitions', nsmap))
        numberOfAcquisitions.text = meta['prod']['numberOfAcquisitions']
    else:
        antennaLookDirection = etree.SubElement(acquisition, _nsc('sar:antennaLookDirection', nsmap))
        antennaLookDirection.text = meta['common']['antennaLookDirection']
        minimumIncidenceAngle = etree.SubElement(acquisition, _nsc('sar:minimumIncidenceAngle', nsmap),
                                                 attrib={'uom': 'deg'})
        minimumIncidenceAngle.text = str(meta['source'][uid]['incidenceAngleMin'])
        maximumIncidenceAngle = etree.SubElement(acquisition, _nsc('sar:maximumIncidenceAngle', nsmap),
                                                 attrib={'uom': 'deg'})
        maximumIncidenceAngle.text = str(meta['source'][uid]['incidenceAngleMax'])
        orbitMeanAltitude = etree.SubElement(acquisition, _nsc('nrb:orbitMeanAltitude', nsmap),
                                             attrib={'uom': 'm'})
        orbitMeanAltitude.text = meta['common']['orbitMeanAltitude']
        dataTakeID = etree.SubElement(acquisition, _nsc('nrb:dataTakeID', nsmap))
        dataTakeID.text = meta['source'][uid]['datatakeID']
        majorCycleID = etree.SubElement(acquisition, _nsc('nrb:majorCycleID', nsmap))
        majorCycleID.text = meta['source'][uid]['majorCycleID']


def om_feature_of_interest(root, nsmap, scene_id, extent, center):
    """
    Creates the `om:featureOfInterest` XML elements.
    
    Parameters
    ----------
    root: etree.Element
        Root XML element.
    nsmap: dict
        Dictionary listing abbreviation (key) and URI (value) of all necessary XML namespaces.
    scene_id: str
        Scene basename.
    extent: str
        Footprint coordinates of the scene.
    center: str
        Center coordinates of the footprint.
    
    Returns
    -------
    None
    """
    featureOfInterest = etree.SubElement(root, _nsc('om:featureOfInterest', nsmap))
    footprint = etree.SubElement(featureOfInterest, _nsc('eop:Footprint', nsmap),
                                 attrib={_nsc('gml:id', nsmap): scene_id + '_5'})
    
    multiExtentOf = etree.SubElement(footprint, _nsc('eop:multiExtentOf', nsmap))
    multiSurface = etree.SubElement(multiExtentOf, _nsc('gml:MultiSurface', nsmap),
                                    attrib={_nsc('gml:id', nsmap): scene_id + '_6'})
    surfaceMember = etree.SubElement(multiSurface, _nsc('gml:surfaceMember', nsmap))
    polygon = etree.SubElement(surfaceMember, _nsc('gml:Polygon', nsmap),
                               attrib={_nsc('gml:id', nsmap): scene_id + '_7'})
    exterior = etree.SubElement(polygon, _nsc('gml:exterior', nsmap))
    linearRing = etree.SubElement(exterior, _nsc('gml:LinearRing', nsmap))
    posList = etree.SubElement(linearRing, _nsc('gml:posList', nsmap))
    posList.text = extent
    
    centerOf = etree.SubElement(footprint, _nsc('eop:centerOf', nsmap))
    point = etree.SubElement(centerOf, _nsc('gml:Point', nsmap), attrib={_nsc('gml:id', nsmap): scene_id + '_8'})
    pos = etree.SubElement(point, _nsc('gml:pos', nsmap))
    pos.text = center


def product_xml(meta, target, tifs, nsmap, exist_ok=False):
    """
    Function to generate product-level metadata for an NRB target product in OGC 10-157r4 compliant XML format.
    
    Parameters
    ----------
    meta: dict
        Metadata dictionary generated with `metadata.extract.meta_dict`
    target: str
        A path pointing to the root directory of a product scene.
    tifs: list[str]
        List of paths to all GeoTIFF files of the currently processed NRB product.
    nsmap: dict
        Dictionary listing abbreviation (key) and URI (value) of all necessary XML namespaces.
    exist_ok: bool
        do not create files if they already exist?
    
    Returns
    -------
    None
    """
    scene_id = os.path.basename(target)
    outname = os.path.join(target, '{}.xml'.format(scene_id))
    if os.path.isfile(outname) and exist_ok:
        return
    print(outname)
    timeCreated = datetime.strftime(meta['prod']['timeCreated'], '%Y-%m-%dT%H:%M:%S.%f')
    timeStart = datetime.strftime(meta['prod']['timeStart'], '%Y-%m-%dT%H:%M:%S.%f')
    timeStop = datetime.strftime(meta['prod']['timeStop'], '%Y-%m-%dT%H:%M:%S.%f')
    
    root = etree.Element(_nsc('nrb:EarthObservation', nsmap), nsmap=nsmap,
                         attrib={_nsc('gml:id', nsmap): scene_id + '_1'})
    om_time(root=root, nsmap=nsmap, scene_id=scene_id, time_start=timeStart, time_stop=timeStop)
    om_procedure(root=root, nsmap=nsmap, scene_id=scene_id, meta=meta, prod=True)
    observedProperty = etree.SubElement(root, _nsc('om:observedProperty', nsmap),
                                        attrib={'nilReason': 'inapplicable'})
    om_feature_of_interest(root=root, nsmap=nsmap, scene_id=scene_id,
                           extent=meta['prod']['geom_xml_envelope'],
                           center=meta['prod']['geom_xml_center'])
    
    ####################################################################################################################
    result = etree.SubElement(root, _nsc('om:result', nsmap))
    earthObservationResult = etree.SubElement(result, _nsc('eop:EarthObservationResult', nsmap),
                                              attrib={_nsc('gml:id', nsmap): scene_id + '_9'})
    product = etree.SubElement(earthObservationResult, _nsc('eop:product', nsmap))
    productInformation = etree.SubElement(product, _nsc('nrb:ProductInformation', nsmap))
    fileName = etree.SubElement(productInformation, _nsc('eop:fileName', nsmap))
    serviceReference = etree.SubElement(fileName, _nsc('ows:ServiceReference', nsmap),
                                        attrib={_nsc('xlink:href', nsmap): scene_id})
    requestMessage = etree.SubElement(serviceReference, _nsc('ows:RequestMessage', nsmap))
    
    for tif in tifs:
        relpath = './' + os.path.relpath(tif, target).replace('\\', '/')
        z_errors = meta['prod']['compression_zerrors']
        pattern = '|'.join(z_errors.keys())
        match = re.search(pattern, os.path.basename(tif))
        
        product = etree.SubElement(earthObservationResult, _nsc('eop:product', nsmap))
        productInformation = etree.SubElement(product, _nsc('nrb:ProductInformation', nsmap))
        fileName = etree.SubElement(productInformation, _nsc('eop:fileName', nsmap))
        serviceReference = etree.SubElement(fileName, _nsc('ows:ServiceReference', nsmap),
                                            attrib={_nsc('xlink:href', nsmap): relpath})
        requestMessage = etree.SubElement(serviceReference, _nsc('ows:RequestMessage', nsmap))
        
        size = etree.SubElement(productInformation, _nsc('eop:size', nsmap), attrib={'uom': 'bytes'})
        size.text = str(os.path.getsize(tif))
        headerSize = etree.SubElement(productInformation, _nsc('nrb:headerSize', nsmap), attrib={'uom': 'bytes'})
        headerSize.text = str(get_header_size(tif))
        byteOrder = etree.SubElement(productInformation, _nsc('nrb:byteOrder', nsmap))
        byteOrder.text = meta['prod']['fileByteOrder']
        dataFormat = etree.SubElement(productInformation, _nsc('nrb:dataFormat', nsmap))
        dataFormat.text = meta['prod']['fileFormat']
        dataType = etree.SubElement(productInformation, _nsc('nrb:dataType', nsmap))
        dataType.text = meta['prod']['fileDataType'].upper()
        bitsPerSample = etree.SubElement(productInformation, _nsc('nrb:bitsPerSample', nsmap))
        bitsPerSample.text = meta['prod']['fileBitsPerSample']
        noDataVal = etree.SubElement(productInformation, _nsc('nrb:noDataValue', nsmap))
        noDataVal.text = 'NaN'
        compressionType = etree.SubElement(productInformation, _nsc('nrb:compressionType', nsmap))
        compressionType.text = meta['prod']['compression_type']
        if match is not None:
            k = match.group()
            compressionzError = etree.SubElement(productInformation, _nsc('nrb:compressionZError', nsmap))
            compressionzError.text = str(z_errors[k])
        
        if 'annotation' in tif:
            key = re.search('-[a-z]{2}(?:-[a-z]{2}|).tif', tif).group()
            np_pat = '-np-[vh]{2}.tif'
            if re.search(np_pat, key) is not None:
                key = np_pat
            
            if key in ['-dm.tif', '-id.tif']:
                dataType.text = 'UINT'
                bitsPerSample.text = '8'
                noDataVal.text = '255'
                
                if key == '-dm.tif':
                    with Raster(tif) as dm_ras:
                        band_descr = [dm_ras.raster.GetRasterBand(band).GetDescription() for band in
                                      range(1, dm_ras.bands + 1)]
                    if 1 < len(band_descr) < len(SAMPLE_MAP[key]['values']):
                        samples = {key: val for key, val in SAMPLE_MAP[key]['values'].items() if val in band_descr}
                        for i, sample_val in enumerate(samples.values()):
                            bitValue = etree.SubElement(productInformation, _nsc('nrb:bitValue', nsmap),
                                                        attrib={'band': str(i + 1),
                                                                'name': sample_val})
                            bitValue.text = '1'
                    else:
                        raise RuntimeError('{} contains an unexpected number of bands!'.format(tif))
                else:  # key == '-id.tif'
                    src_list = list(meta['source'].keys())
                    src_target = [os.path.basename(meta['source'][src]['filename']).replace('.SAFE',
                                                                                            '').replace('.zip', '')
                                  for src in src_list]
                    for i, s in enumerate(src_target):
                        bitValue = etree.SubElement(productInformation, _nsc('nrb:bitValue', nsmap),
                                                    attrib={'band': '1', 'name': s})
                        bitValue.text = str(i + 1)
            
            if SAMPLE_MAP[key]['unit'] is None:
                SAMPLE_MAP[key]['unit'] = 'unitless'
            sampleType = etree.SubElement(productInformation, _nsc('nrb:sampleType', nsmap),
                                          attrib={'uom': SAMPLE_MAP[key]['unit']})
            sampleType.text = SAMPLE_MAP[key]['type']
            
            if key == '-ei.tif':
                ellipsoidalHeight = etree.SubElement(productInformation, _nsc('nrb:ellipsoidalHeight', nsmap),
                                                     attrib={'uom': 'm'})
                ellipsoidalHeight.text = meta['prod']['ellipsoidalHeight']
        
        if 'measurement' in tif:
            creationTime = etree.SubElement(productInformation, _nsc('nrb:creationTime', nsmap))
            creationTime.text = datetime.fromtimestamp(os.path.getctime(tif)).isoformat()
            polarization = etree.SubElement(productInformation, _nsc('nrb:polarization', nsmap))
            polarization.text = re.search('[vh]{2}', tif).group().upper()
            numBorderPixels = etree.SubElement(productInformation, _nsc('nrb:numBorderPixels', nsmap))
            numBorderPixels.text = str(meta['prod']['numBorderPixels'])
    
    ####################################################################################################################
    metaDataProperty = etree.SubElement(root, _nsc('eop:metaDataProperty', nsmap))
    earthObservationMetaData = etree.SubElement(metaDataProperty, _nsc('nrb:EarthObservationMetaData', nsmap))
    
    identifier = etree.SubElement(earthObservationMetaData, _nsc('eop:identifier', nsmap))
    identifier.text = scene_id
    doi = etree.SubElement(earthObservationMetaData, _nsc('eop:doi', nsmap))
    doi.text = meta['prod']['doi']
    acquisitionType = etree.SubElement(earthObservationMetaData, _nsc('eop:acquisitionType', nsmap))
    acquisitionType.text = meta['prod']['acquisitionType']
    status = etree.SubElement(earthObservationMetaData, _nsc('eop:status', nsmap))
    status.text = meta['prod']['status']
    
    processing = etree.SubElement(earthObservationMetaData, _nsc('eop:processing', nsmap))
    processingInformation = etree.SubElement(processing, _nsc('nrb:ProcessingInformation', nsmap))
    processingCenter = etree.SubElement(processingInformation, _nsc('eop:processingCenter', nsmap),
                                        attrib={'codeSpace': 'urn:esa:eop:Sentinel1:facility'})
    processingCenter.text = meta['prod']['processingCenter']
    processingDate = etree.SubElement(processingInformation, _nsc('eop:processingDate', nsmap))
    processingDate.text = timeCreated
    processorName = etree.SubElement(processingInformation, _nsc('eop:processorName', nsmap))
    processorName.text = meta['prod']['processorName']
    processorVersion = etree.SubElement(processingInformation, _nsc('eop:processorVersion', nsmap))
    processorVersion.text = meta['prod']['processorVersion']
    processingMode = etree.SubElement(processingInformation, _nsc('eop:processingMode', nsmap),
                                      attrib={'codeSpace': 'urn:esa:eop:Sentinel1:class'})
    processingMode.text = meta['prod']['processingMode']
    processingLevel = etree.SubElement(processingInformation, _nsc('nrb:processingLevel', nsmap))
    processingLevel.text = meta['common']['processingLevel']
    for src in list(meta['source'].keys()):
        src_path = '{}.xml'.format(os.path.basename(meta['source'][src]['filename']).split('.')[0])
        src_target = os.path.join('./source', src_path).replace('\\', '/')
        sourceProduct = etree.SubElement(processingInformation, _nsc('nrb:sourceProduct', nsmap),
                                         attrib={_nsc('xlink:href', nsmap): src_target})
    auxData1 = etree.SubElement(processingInformation, _nsc('nrb:auxiliaryDataSetFileName', nsmap),
                                attrib={_nsc('xlink:href', nsmap): meta['prod']['ancillaryData_KML']})
    speckleFilterApplied = etree.SubElement(processingInformation, _nsc('nrb:speckleFilterApplied', nsmap))
    speckleFilterApplied.text = str(meta['prod']['speckleFilterApplied']).lower()
    nrApplied = etree.SubElement(processingInformation, _nsc('nrb:NRApplied', nsmap))
    nrApplied.text = str(meta['prod']['NRApplied']).lower()
    if meta['prod']['NRApplied']:
        nrAlgorithm = etree.SubElement(processingInformation, _nsc('nrb:NRAlgorithm', nsmap),
                                       attrib={_nsc('xlink:href', nsmap): meta['prod']['NRAlgorithm']})
    rtcAlgorithm = etree.SubElement(processingInformation, _nsc('nrb:RTCAlgorithm', nsmap),
                                    attrib={_nsc('xlink:href', nsmap): meta['prod']['RTCAlgorithm']})
    geoCorrAlgorithm = etree.SubElement(processingInformation, _nsc('nrb:geoCorrAlgorithm', nsmap),
                                        attrib={_nsc('xlink:href', nsmap): meta['prod']['geoCorrAlgorithm']})
    geoCorrResamplingMethod = etree.SubElement(processingInformation, _nsc('nrb:geoCorrResamplingAlgorithm', nsmap))
    geoCorrResamplingMethod.text = meta['prod']['geoCorrResamplingMethod'].upper()
    demReference = etree.SubElement(processingInformation, _nsc('nrb:DEMReference', nsmap),
                                    attrib={'name': meta['prod']['demName'],
                                            'dem': meta['prod']['demType'],
                                            _nsc('xlink:href', nsmap): meta['prod']['demReference']})
    demResamplingMethod = etree.SubElement(processingInformation, _nsc('nrb:DEMResamplingMethod', nsmap))
    demResamplingMethod.text = meta['prod']['demResamplingMethod'].upper()
    demAccess = etree.SubElement(processingInformation, _nsc('nrb:DEMAccess', nsmap),
                                 attrib={_nsc('xlink:href', nsmap): meta['prod']['demAccess']})
    egmReference = etree.SubElement(processingInformation, _nsc('nrb:EGMReference', nsmap),
                                    attrib={_nsc('xlink:href', nsmap): meta['prod']['demEGMReference']})
    egmResamplingMethod = etree.SubElement(processingInformation, _nsc('nrb:EGMResamplingMethod', nsmap))
    egmResamplingMethod.text = meta['prod']['demEGMResamplingMethod'].upper()
    
    productType = etree.SubElement(earthObservationMetaData, _nsc('nrb:productType', nsmap),
                                   attrib={'codeSpace': 'urn:esa:eop:Sentinel1:class'})
    productType.text = meta['prod']['productName-short']
    azimuthNumberOfLooks = etree.SubElement(earthObservationMetaData, _nsc('nrb:azimuthNumberOfLooks', nsmap))
    azimuthNumberOfLooks.text = str(meta['prod']['azimuthNumberOfLooks'])
    rangeNumberOfLooks = etree.SubElement(earthObservationMetaData, _nsc('nrb:rangeNumberOfLooks', nsmap))
    rangeNumberOfLooks.text = str(meta['prod']['rangeNumberOfLooks'])
    refDoc = etree.SubElement(earthObservationMetaData, _nsc('nrb:refDoc', nsmap),
                              attrib={'name': meta['prod']['productName'],
                                      'version': meta['prod']['card4l-version'],
                                      _nsc('xlink:href', nsmap): meta['prod']['card4l-link']})
    radiometricAccuracyRelative = etree.SubElement(earthObservationMetaData,
                                                   _nsc('nrb:radiometricAccuracyRelative', nsmap), attrib={'uom': 'dB'})
    radiometricAccuracyRelative.text = meta['prod']['radiometricAccuracyRelative']
    radiometricAccuracyAbsolute = etree.SubElement(earthObservationMetaData,
                                                   _nsc('nrb:radiometricAccuracyAbsolute', nsmap), attrib={'uom': 'dB'})
    radiometricAccuracyAbsolute.text = meta['prod']['radiometricAccuracyAbsolute']
    radacc_ref = str(meta['prod']['radiometricAccuracyReference'])
    radiometricAccuracyReference = etree.SubElement(earthObservationMetaData,
                                                    _nsc('nrb:radiometricAccuracyReference', nsmap),
                                                    attrib={_nsc('xlink:href', nsmap): radacc_ref})
    geoCorrAccuracyType = etree.SubElement(earthObservationMetaData, _nsc('nrb:geoCorrAccuracyType', nsmap))
    geoCorrAccuracyType.text = meta['prod']['geoCorrAccuracyType']
    geoCorrAccuracyNorthernSTDev = etree.SubElement(earthObservationMetaData,
                                                    _nsc('nrb:geoCorrAccuracyNorthernSTDev', nsmap),
                                                    attrib={'uom': 'm'})
    geoCorrAccuracyNorthernSTDev.text = meta['prod']['geoCorrAccuracyNorthernSTDev']
    geoCorrAccuracyEasternSTDev = etree.SubElement(earthObservationMetaData,
                                                   _nsc('nrb:geoCorrAccuracyEasternSTDev', nsmap), attrib={'uom': 'm'})
    geoCorrAccuracyEasternSTDev.text = meta['prod']['geoCorrAccuracyEasternSTDev']
    geoCorrAccuracyNorthernBias = etree.SubElement(earthObservationMetaData,
                                                   _nsc('nrb:geoCorrAccuracyNorthernBias', nsmap), attrib={'uom': 'm'})
    geoCorrAccuracyNorthernBias.text = meta['prod']['geoCorrAccuracyNorthernBias']
    geoCorrAccuracyEasternBias = etree.SubElement(earthObservationMetaData,
                                                  _nsc('nrb:geoCorrAccuracyEasternBias', nsmap), attrib={'uom': 'm'})
    geoCorrAccuracyEasternBias.text = meta['prod']['geoCorrAccuracyEasternBias']
    geoCorrAccuracy_rRMSE = etree.SubElement(earthObservationMetaData,
                                             _nsc('nrb:geoCorrAccuracy_rRMSE', nsmap), attrib={'uom': 'm'})
    geoCorrAccuracy_rRMSE.text = meta['prod']['geoCorrAccuracy_rRMSE']
    geoCorrAccuracyReference = etree.SubElement(earthObservationMetaData, _nsc('nrb:geoCorrAccuracyReference', nsmap),
                                                attrib={_nsc('xlink:href', nsmap): meta['prod'][
                                                    'geoCorrAccuracyReference']})
    numLines = etree.SubElement(earthObservationMetaData, _nsc('nrb:numLines', nsmap))
    numLines.text = meta['prod']['numLines']
    numPixelsPerLine = etree.SubElement(earthObservationMetaData, _nsc('nrb:numPixelsPerLine', nsmap))
    numPixelsPerLine.text = meta['prod']['numPixelsPerLine']
    columnSpacing = etree.SubElement(earthObservationMetaData, _nsc('nrb:columnSpacing', nsmap), attrib={'uom': 'm'})
    columnSpacing.text = meta['prod']['pxSpacingColumn']
    rowSpacing = etree.SubElement(earthObservationMetaData, _nsc('nrb:rowSpacing', nsmap), attrib={'uom': 'm'})
    rowSpacing.text = meta['prod']['pxSpacingRow']
    pixelCoordinateConvention = etree.SubElement(earthObservationMetaData, _nsc('nrb:pixelCoordinateConvention', nsmap))
    pixelCoordinateConvention.text = meta['prod']['pixelCoordinateConvention']
    backscatterMeasurement = etree.SubElement(earthObservationMetaData, _nsc('nrb:backscatterMeasurement', nsmap))
    backscatterMeasurement.text = meta['prod']['backscatterMeasurement']
    backscatterConvention = etree.SubElement(earthObservationMetaData, _nsc('nrb:backscatterConvention', nsmap))
    backscatterConvention.text = meta['prod']['backscatterConvention']
    backscatterConversionEq = etree.SubElement(earthObservationMetaData, _nsc('nrb:backscatterConversionEq', nsmap),
                                               attrib={'uom': 'dB'})
    backscatterConversionEq.text = meta['prod']['backscatterConversionEq']
    griddingConvention = etree.SubElement(earthObservationMetaData, _nsc('nrb:griddingConvention', nsmap),
                                          attrib={_nsc('xlink:href', nsmap): meta['prod']['griddingConventionURL']})
    mgrsID = etree.SubElement(earthObservationMetaData, _nsc('nrb:mgrsID', nsmap))
    mgrsID.text = meta['prod']['mgrsID']
    crsEPSG = etree.SubElement(earthObservationMetaData, _nsc('nrb:crsEPSG', nsmap),
                               attrib={'codeSpace': 'urn:esa:eop:crs'})
    crsEPSG.text = meta['prod']['crsEPSG']
    crsWKT = etree.SubElement(earthObservationMetaData, _nsc('nrb:crsWKT', nsmap))
    crsWKT.text = meta['prod']['crsWKT']
    
    ####################################################################################################################
    etree.indent(root)
    tree = etree.ElementTree(root)
    tree.write(outname, pretty_print=True, xml_declaration=True, encoding='utf-8')


def source_xml(meta, target, nsmap, exist_ok=False):
    """
    Function to generate source-level metadata for an NRB target product in OGC 10-157r4 compliant XML format.
    
    Parameters
    ----------
    meta: dict
        Metadata dictionary generated with `metadata.extract.meta_dict`
    target: str
        A path pointing to the root directory of a product scene.
    nsmap: dict
        Dictionary listing abbreviation (key) and URI (value) of all necessary XML namespaces.
    exist_ok: bool
        do not create files if they already exist?
    
    Returns
    -------
    None
    """
    metadir = os.path.join(target, 'source')
    os.makedirs(metadir, exist_ok=True)
    
    for uid in list(meta['source'].keys()):
        scene = os.path.basename(meta['source'][uid]['filename']).split('.')[0]
        outname = os.path.join(metadir, '{}.xml'.format(scene))
        if os.path.isfile(outname) and exist_ok:
            continue
        print(outname)
        timeStart = datetime.strftime(meta['source'][uid]['timeStart'], '%Y-%m-%dT%H:%M:%S.%f')
        timeStop = datetime.strftime(meta['source'][uid]['timeStop'], '%Y-%m-%dT%H:%M:%S.%f')
        
        root = etree.Element(_nsc('nrb:EarthObservation', nsmap), nsmap=nsmap,
                             attrib={_nsc('gml:id', nsmap): scene + '_1'})
        om_time(root=root, nsmap=nsmap, scene_id=scene, time_start=timeStart, time_stop=timeStop)
        om_procedure(root=root, nsmap=nsmap, scene_id=scene, meta=meta, uid=uid, prod=False)
        observedProperty = etree.SubElement(root, _nsc('om:observedProperty', nsmap),
                                            attrib={'nilReason': 'inapplicable'})
        om_feature_of_interest(root=root, nsmap=nsmap, scene_id=scene,
                               extent=meta['source'][uid]['geom_xml_envelop'],
                               center=meta['source'][uid]['geom_xml_center'])
        
        ################################################################################################################
        result = etree.SubElement(root, _nsc('om:result', nsmap))
        earthObservationResult = etree.SubElement(result, _nsc('eop:EarthObservationResult', nsmap),
                                                  attrib={_nsc('gml:id', nsmap): scene + '_9'})
        product = etree.SubElement(earthObservationResult, _nsc('eop:product', nsmap))
        productInformation = etree.SubElement(product, _nsc('nrb:ProductInformation', nsmap))
        fileName = etree.SubElement(productInformation, _nsc('eop:fileName', nsmap))
        serviceReference = etree.SubElement(fileName, _nsc('ows:ServiceReference', nsmap),
                                            attrib={_nsc('xlink:href', nsmap): scene})
        requestMessage = etree.SubElement(serviceReference, _nsc('ows:RequestMessage', nsmap))
        
        ################################################################################################################
        metaDataProperty = etree.SubElement(root, _nsc('eop:metaDataProperty', nsmap))
        earthObservationMetaData = etree.SubElement(metaDataProperty, _nsc('nrb:EarthObservationMetaData', nsmap))
        
        identifier = etree.SubElement(earthObservationMetaData, _nsc('eop:identifier', nsmap))
        identifier.text = scene
        doi = etree.SubElement(earthObservationMetaData, _nsc('eop:doi', nsmap))
        doi.text = meta['source'][uid]['doi']
        acquisitionType = etree.SubElement(earthObservationMetaData, _nsc('eop:acquisitionType', nsmap))
        acquisitionType.text = meta['source'][uid]['acquisitionType']
        status = etree.SubElement(earthObservationMetaData, _nsc('eop:status', nsmap))
        status.text = meta['source'][uid]['status']
        
        processing = etree.SubElement(earthObservationMetaData, _nsc('eop:processing', nsmap))
        processingInformation = etree.SubElement(processing, _nsc('nrb:ProcessingInformation', nsmap))
        processingCenter = etree.SubElement(processingInformation, _nsc('eop:processingCenter', nsmap),
                                            attrib={'codeSpace': 'urn:esa:eop:Sentinel1:facility'})
        processingCenter.text = meta['source'][uid]['processingCenter']
        processingDate = etree.SubElement(processingInformation, _nsc('eop:processingDate', nsmap))
        processingDate.text = meta['source'][uid]['processingDate']
        processorName = etree.SubElement(processingInformation, _nsc('eop:processorName', nsmap))
        processorName.text = meta['source'][uid]['processorName']
        processorVersion = etree.SubElement(processingInformation, _nsc('eop:processorVersion', nsmap))
        processorVersion.text = meta['source'][uid]['processorVersion']
        processingMode = etree.SubElement(processingInformation, _nsc('eop:processingMode', nsmap))
        processingMode.text = meta['source'][uid]['processingMode']
        processingLevel = etree.SubElement(processingInformation, _nsc('nrb:processingLevel', nsmap))
        processingLevel.text = meta['common']['processingLevel']
        orbitDataSource = etree.SubElement(processingInformation, _nsc('nrb:orbitDataSource', nsmap))
        orbitDataSource.text = meta['source'][uid]['orbitDataSource'].upper()
        orbitStateVector = etree.SubElement(processingInformation, _nsc('nrb:orbitStateVector', nsmap),
                                            attrib={'access': meta['source'][uid]['orbitDataAccess']})
        orbitStateVector.text = meta['source'][uid]['orbitStateVector']
        for swath in meta['source'][uid]['swaths']:
            azimuthLookBandwidth = etree.SubElement(processingInformation, _nsc('nrb:azimuthLookBandwidth', nsmap),
                                                    attrib={'uom': 'Hz', 'beam': swath})
            azimuthLookBandwidth.text = str(meta['source'][uid]['azimuthLookBandwidth'][swath])
        for swath in meta['source'][uid]['swaths']:
            rangeLookBandwidth = etree.SubElement(processingInformation, _nsc('nrb:rangeLookBandwidth', nsmap),
                                                  attrib={'uom': 'Hz', 'beam': swath})
            rangeLookBandwidth.text = str(meta['source'][uid]['rangeLookBandwidth'][swath])
        lutApplied = etree.SubElement(processingInformation, _nsc('nrb:lutApplied', nsmap))
        lutApplied.text = meta['source'][uid]['lutApplied']
        
        productType = etree.SubElement(earthObservationMetaData, _nsc('nrb:productType', nsmap),
                                       attrib={'codeSpace': 'urn:esa:eop:Sentinel1:class'})
        productType.text = meta['source'][uid]['productType']
        for swath in meta['source'][uid]['swaths']:
            azimuthNumberOfLooks = etree.SubElement(earthObservationMetaData,
                                                    _nsc('nrb:azimuthNumberOfLooks', nsmap),
                                                    attrib={'beam': swath})
            azimuthNumberOfLooks.text = meta['source'][uid]['azimuthNumberOfLooks'][swath]
        for swath in meta['source'][uid]['swaths']:
            rangeNumberOfLooks = etree.SubElement(earthObservationMetaData,
                                                  _nsc('nrb:rangeNumberOfLooks', nsmap),
                                                  attrib={'beam': swath})
            rangeNumberOfLooks.text = meta['source'][uid]['rangeNumberOfLooks'][swath]
        dataGeometry = etree.SubElement(earthObservationMetaData,
                                        _nsc('nrb:dataGeometry', nsmap))
        dataGeometry.text = meta['source'][uid]['dataGeometry']
        for swath in meta['source'][uid]['swaths']:
            azimuthResolution = etree.SubElement(earthObservationMetaData,
                                                 _nsc('nrb:azimuthResolution', nsmap),
                                                 attrib={'uom': 'm', 'beam': swath})
            azimuthResolution.text = str(meta['source'][uid]['azimuthResolution'][swath])
        for swath in meta['source'][uid]['swaths']:
            rangeResolution = etree.SubElement(earthObservationMetaData,
                                               _nsc('nrb:rangeResolution', nsmap),
                                               attrib={'uom': 'm', 'beam': swath})
            rangeResolution.text = str(meta['source'][uid]['rangeResolution'][swath])
        azimuthPixelSpacing = etree.SubElement(earthObservationMetaData, _nsc('nrb:azimuthPixelSpacing', nsmap),
                                               attrib={'uom': 'm'})
        azimuthPixelSpacing.text = str(mean(meta['source'][uid]['azimuthPixelSpacing'].values()))
        rangePixelSpacing = etree.SubElement(earthObservationMetaData, _nsc('nrb:rangePixelSpacing', nsmap),
                                             attrib={'uom': 'm'})
        rangePixelSpacing.text = str(mean(meta['source'][uid]['rangePixelSpacing'].values()))
        
        performance = etree.SubElement(earthObservationMetaData, _nsc('nrb:performance', nsmap))
        performanceIndicators = etree.SubElement(performance, _nsc('nrb:PerformanceIndicators', nsmap))
        noiseEquivalentIntensityType = etree.SubElement(performanceIndicators,
                                                        _nsc('nrb:noiseEquivalentIntensityType', nsmap),
                                                        attrib={'uom': 'dB'})
        noiseEquivalentIntensityType.text = str(meta['source'][uid]['perfNoiseEquivalentIntensityType'])
        for pol in meta['common']['polarisationChannels']:
            estimatesMin = etree.SubElement(performanceIndicators, _nsc('nrb:estimates', nsmap),
                                            attrib={'pol': pol, 'type': 'minimum'})
            estimatesMin.text = str(meta['source'][uid]['perfEstimates'][pol]['minimum'])
            estimatesMax = etree.SubElement(performanceIndicators, _nsc('nrb:estimates', nsmap),
                                            attrib={'pol': pol, 'type': 'maximum'})
            estimatesMax.text = str(meta['source'][uid]['perfEstimates'][pol]['maximum'])
            estimatesMean = etree.SubElement(performanceIndicators, _nsc('nrb:estimates', nsmap),
                                             attrib={'pol': pol, 'type': 'mean'})
            estimatesMean.text = str(meta['source'][uid]['perfEstimates'][pol]['mean'])
        equivalentNumberOfLooks = etree.SubElement(performanceIndicators, _nsc('nrb:equivalentNumberOfLooks', nsmap))
        equivalentNumberOfLooks.text = str(meta['source'][uid]['perfEquivalentNumberOfLooks'])
        peakSideLobeRatio = etree.SubElement(performanceIndicators, _nsc('nrb:peakSideLobeRatio', nsmap),
                                             attrib={'uom': 'dB'})
        peakSideLobeRatio.text = str(meta['source'][uid]['perfPeakSideLobeRatio'])
        integratedSideLobeRatio = etree.SubElement(performanceIndicators, _nsc('nrb:integratedSideLobeRatio', nsmap),
                                                   attrib={'uom': 'dB'})
        integratedSideLobeRatio.text = str(meta['source'][uid]['perfIntegratedSideLobeRatio'])
        
        polCalMatrices = etree.SubElement(earthObservationMetaData, _nsc('nrb:polCalMatrices', nsmap),
                                          attrib={
                                              _nsc('xlink:href', nsmap): str(meta['source'][uid]['polCalMatrices'])})
        meanFaradayRotationAngle = etree.SubElement(earthObservationMetaData,
                                                    _nsc('nrb:meanFaradayRotationAngle', nsmap), attrib={'uom': 'deg'})
        meanFaradayRotationAngle.text = meta['source'][uid]['faradayMeanRotationAngle']
        faraday_ref = str(meta['source'][uid]['faradayRotationReference'])
        referenceFaradayRotation = etree.SubElement(earthObservationMetaData,
                                                    _nsc('nrb:referenceFaradayRotation', nsmap),
                                                    attrib={_nsc('xlink:href', nsmap): faraday_ref})
        ionosphereIndicator = etree.SubElement(earthObservationMetaData, _nsc('nrb:ionosphereIndicator', nsmap))
        ionosphereIndicator.text = meta['source'][uid]['ionosphereIndicator']
        
        ################################################################################################################
        etree.indent(root)
        tree = etree.ElementTree(root)
        tree.write(outname, pretty_print=True, xml_declaration=True, encoding='utf-8')


def main(meta, target, tifs, exist_ok=False):
    """
    Wrapper for `source_xml` and `product_xml`.
    
    Parameters
    ----------
    meta: dict
        Metadata dictionary generated with `metadata.extract.meta_dict`
    target: str
        A path pointing to the root directory of a product scene.
    tifs: list[str]
        List of paths to all GeoTIFF files of the currently processed NRB product.
    exist_ok: bool
        do not create files if they already exist?
    
    Returns
    -------
    None
    """
    NS_MAP_prod = deepcopy(NS_MAP)
    NS_MAP_src = deepcopy(NS_MAP)
    NS_MAP_prod['nrb'] = NS_MAP['nrb']['product']
    NS_MAP_src['nrb'] = NS_MAP['nrb']['source']
    
    source_xml(meta=meta, target=target, nsmap=NS_MAP_src, exist_ok=exist_ok)
    product_xml(meta=meta, target=target, tifs=tifs, nsmap=NS_MAP_prod, exist_ok=exist_ok)
