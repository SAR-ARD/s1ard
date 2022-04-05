import os
import re
from lxml import etree
from datetime import datetime
from spatialist import Raster
from S1_NRB.metadata.mapping import SAMPLE_MAP, NS_MAP
from S1_NRB.metadata.extract import get_header_size


def _nsc(text):
    ns, key = text.split(':')
    return '{{{0}}}{1}'.format(NS_MAP[ns], key)


def _common_procedure_elements(eo_equipment, meta, prod=True):
    """
    Adds common (source & product) XML subelements to the `om:procedure/eop:earthObservationEquipment` element.
    
    Parameters
    ----------
    eo_equipment: etree.Element
        `eop:earthObservationEquipment` XML subelement of `om:procedure`, which is one of the main properties of the
        root XML element.
    meta: dict
        Metadata dictionary generated with `metadata.extract.meta_dict`
    prod: bool, optional
        Return XML subelements for further usage in `product_xml` parsing function? Default is True. If False, the
        XML subelements for further usage in the `source_xml` parsing function will be returned.
    
    Returns
    -------
    etree.Element
    """
    earthObservationEquipment = eo_equipment
    
    platform0 = etree.SubElement(earthObservationEquipment, _nsc('eop:platform'))
    platform1 = etree.SubElement(platform0, _nsc('nrb:Platform'))
    shortName = etree.SubElement(platform1, _nsc('eop:shortName'))
    shortName.text = meta['common']['platformShortName'].upper()
    serialIdentifier = etree.SubElement(platform1, _nsc('eop:serialIdentifier'))
    serialIdentifier.text = meta['common']['platformIdentifier']
    instrument0 = etree.SubElement(earthObservationEquipment, _nsc('eop:instrument'))
    instrument1 = etree.SubElement(instrument0, _nsc('eop:Instrument'))
    shortName = etree.SubElement(instrument1, _nsc('eop:shortName'))
    shortName.text = meta['common']['instrumentShortName']
    sensor0 = etree.SubElement(earthObservationEquipment, _nsc('eop:sensor'))
    sensor1 = etree.SubElement(sensor0, _nsc('nrb:Sensor'))
    sensorType = etree.SubElement(sensor1, _nsc('eop:sensorType'))
    sensorType.text = meta['common']['sensorType']
    radarBand = etree.SubElement(sensor1, _nsc('nrb:radarBand'))
    radarBand.text = meta['common']['radarBand']
    operationalMode = etree.SubElement(sensor1, _nsc('eop:operationalMode'),
                                       attrib={'codeSpace': 'urn:esa:eop:C-SAR:operationalMode'})
    operationalMode.text = meta['common']['operationalMode']
    swathIdentifier = etree.SubElement(sensor1, _nsc('eop:swathIdentifier'),
                                       attrib={'codeSpace': 'urn:esa:eop:C-SAR:swathIdentifier'})
    swathIdentifier.text = meta['common']['swathIdentifier']
    acquisitionParameters = etree.SubElement(earthObservationEquipment, _nsc('eop:acquisitionParameters'))
    acquisition = etree.SubElement(acquisitionParameters, _nsc('nrb:Acquisition'))
    polarisationMode = etree.SubElement(acquisition, _nsc('sar:polarisationMode'))
    polarisationMode.text = meta['common']['polarisationMode']
    polarisationChannels = etree.SubElement(acquisition, _nsc('sar:polarisationChannels'))
    polarisationChannels.text = ', '.join(meta['common']['polarisationChannels'])
    orbitDirection = etree.SubElement(acquisition, _nsc('eop:orbitDirection'))
    orbitDirection.text = meta['common']['orbitDirection'].upper()
    orbitNumber = etree.SubElement(acquisition, _nsc('eop:orbitNumber'))
    orbitNumber.text = meta['common']['orbitNumber']
    wrsLongitudeGrid = etree.SubElement(acquisition, _nsc('eop:wrsLongitudeGrid'),
                                        attrib={'codeSpace': 'urn:esa:eop:Sentinel1:relativeOrbits'})
    wrsLongitudeGrid.text = meta['common']['wrsLongitudeGrid']
    
    if prod:
        return acquisition
    else:
        return platform1, sensor1, acquisition


def product_xml(meta, target, tifs):
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
    
    Returns
    -------
    None
    """
    scene_id = os.path.basename(target)
    outname = os.path.join(target, '{}.xml'.format(scene_id))
    print(outname)
    timeCreated = datetime.strftime(meta['prod']['timeCreated'], '%Y-%m-%dT%H:%M:%S.%f')
    timeStart = datetime.strftime(meta['prod']['timeStart'], '%Y-%m-%dT%H:%M:%S.%f')
    timeStop = datetime.strftime(meta['prod']['timeStop'], '%Y-%m-%dT%H:%M:%S.%f')
    
    root = etree.Element(_nsc('nrb:EarthObservation'), nsmap=NS_MAP,
                         attrib={_nsc('gml:id'): scene_id + '_1'})
    
    ####################################################################################################################
    phenomenonTime = etree.SubElement(root, _nsc('om:phenomenonTime'))
    timePeriod = etree.SubElement(phenomenonTime, _nsc('gml:TimePeriod'), attrib={_nsc('gml:id'): scene_id + '_2'})
    beginPosition = etree.SubElement(timePeriod, _nsc('gml:beginPosition'))
    beginPosition.text = timeStart
    endPosition = etree.SubElement(timePeriod, _nsc('gml:endPosition'))
    endPosition.text = timeStop
    
    resultTime = etree.SubElement(root, _nsc('om:resultTime'))
    timeInstant = etree.SubElement(resultTime, _nsc('gml:TimeInstant'), attrib={_nsc('gml:id'): scene_id + '_3'})
    timePosition = etree.SubElement(timeInstant, _nsc('gml:timePosition'))
    timePosition.text = timeStop
    
    ####################################################################################################################
    procedure = etree.SubElement(root, _nsc('om:procedure'))
    earthObservationEquipment = etree.SubElement(procedure, _nsc('eop:EarthObservationEquipment'),
                                                 attrib={_nsc('gml:id'): scene_id + '_4'})
    acquisition = _common_procedure_elements(eo_equipment=earthObservationEquipment, meta=meta, prod=True)
    
    numberOfAcquisitions = etree.SubElement(acquisition, _nsc('nrb:numberOfAcquisitions'))
    numberOfAcquisitions.text = meta['prod']['numberOfAcquisitions']
    
    ####################################################################################################################
    observedProperty = etree.SubElement(root, _nsc('om:observedProperty'), attrib={'nilReason': 'inapplicable'})
    
    ####################################################################################################################
    featureOfInterest = etree.SubElement(root, _nsc('om:featureOfInterest'))
    footprint = etree.SubElement(featureOfInterest, _nsc('eop:Footprint'), attrib={_nsc('gml:id'): scene_id + '_5'})
    
    multiExtentOf = etree.SubElement(footprint, _nsc('eop:multiExtentOf'))
    multiSurface = etree.SubElement(multiExtentOf, _nsc('gml:MultiSurface'), attrib={_nsc('gml:id'): scene_id + '_6'})
    surfaceMember = etree.SubElement(multiSurface, _nsc('gml:surfaceMember'))
    polygon = etree.SubElement(surfaceMember, _nsc('gml:Polygon'), attrib={_nsc('gml:id'): scene_id + '_7'})
    exterior = etree.SubElement(polygon, _nsc('gml:exterior'))
    linearRing = etree.SubElement(exterior, _nsc('gml:LinearRing'))
    posList = etree.SubElement(linearRing, _nsc('gml:posList'), attrib={'uom': 'deg'})
    posList.text = meta['prod']['geom_xml_envelope']
    
    centerOf = etree.SubElement(footprint, _nsc('eop:centerOf'))
    point = etree.SubElement(centerOf, _nsc('gml:Point'), attrib={_nsc('gml:id'): scene_id + '_8'})
    pos = etree.SubElement(point, _nsc('gml:pos'), attrib={'uom': 'deg'})
    pos.text = meta['prod']['geom_xml_center']
    
    ####################################################################################################################
    result = etree.SubElement(root, _nsc('om:result'))
    earthObservationResult = etree.SubElement(result, _nsc('eop:EarthObservationResult'),
                                              attrib={_nsc('gml:id'): scene_id + '_9'})
    product = etree.SubElement(earthObservationResult, _nsc('eop:product'))
    productInformation = etree.SubElement(product, _nsc('nrb:ProductInformation'))
    fileName = etree.SubElement(productInformation, _nsc('eop:fileName'))
    serviceReference = etree.SubElement(fileName, _nsc('ows:ServiceReference'), attrib={_nsc('xlink:href'): scene_id})
    requestMessage = etree.SubElement(serviceReference, _nsc('ows:RequestMessage'))
    
    for tif in tifs:
        relpath = './' + os.path.relpath(tif, target).replace('\\', '/')
        z_errors = meta['prod']['compression_zerrors']
        pattern = '|'.join(z_errors.keys())
        match = re.search(pattern, os.path.basename(tif))
        
        product = etree.SubElement(earthObservationResult, _nsc('eop:product'))
        productInformation = etree.SubElement(product, _nsc('nrb:ProductInformation'))
        fileName = etree.SubElement(productInformation, _nsc('eop:fileName'))
        serviceReference = etree.SubElement(fileName, _nsc('ows:ServiceReference'),
                                            attrib={_nsc('xlink:href'): relpath})
        requestMessage = etree.SubElement(serviceReference, _nsc('ows:RequestMessage'))
        
        size = etree.SubElement(productInformation, _nsc('eop:size'), attrib={'uom': 'bytes'})
        size.text = str(os.path.getsize(tif))
        headerSize = etree.SubElement(productInformation, _nsc('nrb:headerSize'), attrib={'uom': 'bytes'})
        headerSize.text = str(get_header_size(tif))
        byteOrder = etree.SubElement(productInformation, _nsc('nrb:byteOrder'))
        byteOrder.text = meta['prod']['fileByteOrder']
        dataFormat = etree.SubElement(productInformation, _nsc('nrb:dataFormat'))
        dataFormat.text = meta['prod']['fileFormat']
        dataType = etree.SubElement(productInformation, _nsc('nrb:dataType'))
        dataType.text = meta['prod']['fileDataType'].upper()
        bitsPerSample = etree.SubElement(productInformation, _nsc('nrb:bitsPerSample'))
        bitsPerSample.text = meta['prod']['fileBitsPerSample']
        noDataVal = etree.SubElement(productInformation, _nsc('nrb:noDataValue'))
        noDataVal.text = 'NaN'
        compressionType = etree.SubElement(productInformation, _nsc('nrb:compressionType'))
        compressionType.text = meta['prod']['compression_type']
        if match is not None:
            k = match.group()
            compressionzError = etree.SubElement(productInformation, _nsc('nrb:compressionZError'))
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
                            bitValue = etree.SubElement(productInformation, _nsc('nrb:bitValue'),
                                                        attrib={'band': str(i + 1),
                                                                'name': sample_val})
                            bitValue.text = '1'
                    else:
                        raise RuntimeError('{} contains an unexpected number of bands!'.format(tif))
                else:  # key == '-id.tif'
                    src_list = list(meta['source'].keys())
                    src_target = [os.path.basename(meta['source'][src]['filename']).replace('.SAFE', '').replace('.zip', '')
                                  for src in src_list]
                    for i, s in enumerate(src_target):
                        bitValue = etree.SubElement(productInformation, _nsc('nrb:bitValue'),
                                                    attrib={'band': '1', 'name': s})
                        bitValue.text = str(i+1)
            
            if key == '-ei.tif':
                ellipsoidalHeight = etree.SubElement(productInformation, _nsc('nrb:ellipsoidalHeight'),
                                                     attrib={'uom': 'm'})
                ellipsoidalHeight.text = meta['prod']['ellipsoidalHeight']
            
            if SAMPLE_MAP[key]['unit'] is None:
                SAMPLE_MAP[key]['unit'] = 'unitless'
            sampleType = etree.SubElement(productInformation, _nsc('nrb:sampleType'),
                                          attrib={'uom': SAMPLE_MAP[key]['unit']})
            sampleType.text = SAMPLE_MAP[key]['type']
        
        if 'measurement' in tif:
            creationTime = etree.SubElement(productInformation, _nsc('nrb:creationTime'))
            creationTime.text = datetime.fromtimestamp(os.path.getctime(tif)).isoformat()
            polarization = etree.SubElement(productInformation, _nsc('nrb:polarization'))
            polarization.text = re.search('[vh]{2}', tif).group().upper()
            numBorderPixels = etree.SubElement(productInformation, _nsc('nrb:numBorderPixels'))
            numBorderPixels.text = str(meta['prod']['numBorderPixels'])
    
    ####################################################################################################################
    metaDataProperty = etree.SubElement(root, _nsc('eop:metaDataProperty'))
    earthObservationMetaData = etree.SubElement(metaDataProperty, _nsc('nrb:EarthObservationMetaData'))
    
    identifier = etree.SubElement(earthObservationMetaData, _nsc('eop:identifier'))
    identifier.text = scene_id
    doi = etree.SubElement(earthObservationMetaData, _nsc('eop:doi'))
    doi.text = meta['prod']['doi']
    status = etree.SubElement(earthObservationMetaData, _nsc('eop:status'))
    status.text = meta['prod']['status']
    acquisitionType = etree.SubElement(earthObservationMetaData, _nsc('eop:acquisitionType'))
    acquisitionType.text = meta['prod']['acquisitionType']
    productType = etree.SubElement(earthObservationMetaData, _nsc('nrb:productType'),
                                   attrib={'codeSpace': 'urn:esa:eop:Sentinel1:class'})
    productType.text = meta['prod']['productName-short']
    refDoc = etree.SubElement(earthObservationMetaData, _nsc('nrb:refDoc'),
                              attrib={'name': meta['prod']['productName'],
                                      'version': meta['prod']['card4l-version'],
                                      _nsc('xlink:href'): meta['prod']['card4l-link']})
    
    processing = etree.SubElement(earthObservationMetaData, _nsc('eop:processing'))
    processingInformation = etree.SubElement(processing, _nsc('nrb:ProcessingInformation'))
    processingCenter = etree.SubElement(processingInformation, _nsc('eop:processingCenter'),
                                        attrib={'codeSpace': 'urn:esa:eop:Sentinel1:facility'})
    processingCenter.text = meta['prod']['processingCenter']
    processingDate = etree.SubElement(processingInformation, _nsc('eop:processingDate'))
    processingDate.text = timeCreated
    processorName = etree.SubElement(processingInformation, _nsc('eop:processorName'))
    processorName.text = meta['prod']['processorName']
    processorVersion = etree.SubElement(processingInformation, _nsc('eop:processorVersion'))
    processorVersion.text = meta['prod']['processorVersion']
    processingLevel = etree.SubElement(processingInformation, _nsc('eop:processingLevel'))
    processingLevel.text = meta['common']['processingLevel']
    processingMode = etree.SubElement(processingInformation, _nsc('eop:processingMode'),
                                      attrib={'codeSpace': 'urn:esa:eop:Sentinel1:class'})
    processingMode.text = meta['prod']['processingMode']
    for src in list(meta['source'].keys()):
        src_path = '{}.xml'.format(os.path.basename(meta['source'][src]['filename']).split('.')[0])
        src_target = os.path.join('./source', src_path).replace('\\', '/')
        sourceProduct = etree.SubElement(processingInformation, _nsc('nrb:sourceProduct'),
                                         attrib={_nsc('xlink:href'): src_target})
    auxData1 = etree.SubElement(processingInformation, _nsc('nrb:auxiliaryDataSetFileName'),
                                attrib={_nsc('xlink:href'): meta['prod']['ancillaryData_KML']})
    speckleFilterApplied = etree.SubElement(processingInformation, _nsc('nrb:speckleFilterApplied'))
    speckleFilterApplied.text = str(meta['prod']['speckleFilterApplied']).lower()
    nrApplied = etree.SubElement(processingInformation, _nsc('nrb:NRApplied'))
    nrApplied.text = str(meta['prod']['NRApplied']).lower()
    if meta['prod']['NRApplied']:
        nrAlgorithm = etree.SubElement(processingInformation, _nsc('nrb:NRAlgorithm'),
                                       attrib={_nsc('xlink:href'): meta['prod']['NRAlgorithm']})
    rtcAlgorithm = etree.SubElement(processingInformation, _nsc('nrb:RTCAlgorithm'),
                                    attrib={_nsc('xlink:href'): meta['prod']['RTCAlgorithm']})
    geoCorrAlgorithm = etree.SubElement(processingInformation, _nsc('nrb:geoCorrAlgorithm'),
                                        attrib={_nsc('xlink:href'): meta['prod']['geoCorrAlgorithm']})
    geoCorrResamplingMethod = etree.SubElement(processingInformation, _nsc('nrb:geoCorrResamplingAlgorithm'))
    geoCorrResamplingMethod.text = meta['prod']['geoCorrResamplingMethod'].upper()
    demReference = etree.SubElement(processingInformation, _nsc('nrb:DEMReference'),
                                    attrib={'name': meta['prod']['demName'],
                                            'dem': meta['prod']['demType'],
                                            _nsc('xlink:href'): meta['prod']['demReference']})
    demResamplingMethod = etree.SubElement(processingInformation, _nsc('nrb:DEMResamplingMethod'))
    demResamplingMethod.text = meta['prod']['demResamplingMethod'].upper()
    demAccess = etree.SubElement(processingInformation, _nsc('nrb:DEMAccess'),
                                 attrib={_nsc('xlink:href'): meta['prod']['demAccess']})
    egmReference = etree.SubElement(processingInformation, _nsc('nrb:EGMReference'),
                                    attrib={_nsc('xlink:href'): meta['prod']['demEGMReference']})
    egmResamplingMethod = etree.SubElement(processingInformation, _nsc('nrb:EGMResamplingMethod'))
    egmResamplingMethod.text = meta['prod']['demEGMResamplingMethod'].upper()
    
    radiometricAccuracyRelative = etree.SubElement(earthObservationMetaData, _nsc('nrb:radiometricAccuracyRelative'),
                                                   attrib={'uom': 'dB'})
    radiometricAccuracyRelative.text = meta['prod']['radiometricAccuracyRelative']
    radiometricAccuracyAbsolute = etree.SubElement(earthObservationMetaData, _nsc('nrb:radiometricAccuracyAbsolute'),
                                                   attrib={'uom': 'dB'})
    radiometricAccuracyAbsolute.text = meta['prod']['radiometricAccuracyAbsolute']
    radacc_ref = str(meta['prod']['radiometricAccuracyReference'])
    radiometricAccuracyReference = etree.SubElement(earthObservationMetaData, _nsc('nrb:radiometricAccuracyReference'),
                                                    attrib={_nsc('xlink:href'): radacc_ref})
    geoCorrAccuracyType = etree.SubElement(earthObservationMetaData, _nsc('nrb:geoCorrAccuracyType'))
    geoCorrAccuracyType.text = meta['prod']['geoCorrAccuracyType']
    geoCorrAccuracyNorthernSTDev = etree.SubElement(earthObservationMetaData, _nsc('nrb:geoCorrAccuracyNorthernSTDev'),
                                                    attrib={'uom': 'm'})
    geoCorrAccuracyNorthernSTDev.text = meta['prod']['geoCorrAccuracyNorthernSTDev']
    geoCorrAccuracyEasternSTDev = etree.SubElement(earthObservationMetaData, _nsc('nrb:geoCorrAccuracyEasternSTDev'),
                                                   attrib={'uom': 'm'})
    geoCorrAccuracyEasternSTDev.text = meta['prod']['geoCorrAccuracyEasternSTDev']
    geoCorrAccuracyNorthernBias = etree.SubElement(earthObservationMetaData, _nsc('nrb:geoCorrAccuracyNorthernBias'),
                                                   attrib={'uom': 'm'})
    geoCorrAccuracyNorthernBias.text = meta['prod']['geoCorrAccuracyNorthernBias']
    geoCorrAccuracyEasternBias = etree.SubElement(earthObservationMetaData, _nsc('nrb:geoCorrAccuracyEasternBias'),
                                                  attrib={'uom': 'm'})
    geoCorrAccuracyEasternBias.text = meta['prod']['geoCorrAccuracyEasternBias']
    geoCorrAccuracy_rRMSE = etree.SubElement(earthObservationMetaData, _nsc('nrb:geoCorrAccuracy_rRMSE'),
                                             attrib={'uom': 'm'})
    geoCorrAccuracy_rRMSE.text = meta['prod']['geoCorrAccuracy_rRMSE']
    geoCorrAccuracyReference = etree.SubElement(earthObservationMetaData, _nsc('nrb:geoCorrAccuracyReference'),
                                                attrib={_nsc('xlink:href'): meta['prod']['geoCorrAccuracyReference']})
    azimuthNumberOfLooks = etree.SubElement(earthObservationMetaData, _nsc('nrb:azimuthNumberOfLooks'))
    azimuthNumberOfLooks.text = meta['prod']['azimuthNumberOfLooks']
    rangeNumberOfLooks = etree.SubElement(earthObservationMetaData, _nsc('nrb:rangeNumberOfLooks'))
    rangeNumberOfLooks.text = meta['prod']['rangeNumberOfLooks']
    numLines = etree.SubElement(earthObservationMetaData, _nsc('nrb:numLines'))
    numLines.text = meta['prod']['numLines']
    numPixelsPerLine = etree.SubElement(earthObservationMetaData, _nsc('nrb:numPixelsPerLine'))
    numPixelsPerLine.text = meta['prod']['numPixelsPerLine']
    columnSpacing = etree.SubElement(earthObservationMetaData, _nsc('nrb:columnSpacing'), attrib={'uom': 'm'})
    columnSpacing.text = meta['prod']['pxSpacingColumn']
    rowSpacing = etree.SubElement(earthObservationMetaData, _nsc('nrb:rowSpacing'), attrib={'uom': 'm'})
    rowSpacing.text = meta['prod']['pxSpacingRow']
    pixelCoordinateConvention = etree.SubElement(earthObservationMetaData, _nsc('nrb:pixelCoordinateConvention'))
    pixelCoordinateConvention.text = meta['prod']['pixelCoordinateConvention']
    backscatterMeasurement = etree.SubElement(earthObservationMetaData, _nsc('nrb:backscatterMeasurement'))
    backscatterMeasurement.text = meta['prod']['backscatterMeasurement']
    backscatterConvention = etree.SubElement(earthObservationMetaData, _nsc('nrb:backscatterConvention'))
    backscatterConvention.text = meta['prod']['backscatterConvention']
    backscatterConversionEq = etree.SubElement(earthObservationMetaData, _nsc('nrb:backscatterConversionEq'),
                                               attrib={'uom': 'dB'})
    backscatterConversionEq.text = meta['prod']['backscatterConversionEq']
    griddingConvention = etree.SubElement(earthObservationMetaData, _nsc('nrb:griddingConvention'),
                                          attrib={_nsc('xlink:href'): meta['prod']['griddingConventionURL']})
    mgrsID = etree.SubElement(earthObservationMetaData, _nsc('nrb:mgrsID'))
    mgrsID.text = meta['prod']['mgrsID']
    crsEPSG = etree.SubElement(earthObservationMetaData, _nsc('nrb:crsEPSG'), attrib={'codeSpace': 'urn:esa:eop:crs'})
    crsEPSG.text = meta['prod']['crsEPSG']
    crsWKT = etree.SubElement(earthObservationMetaData, _nsc('nrb:crsWKT'))
    crsWKT.text = meta['prod']['crsWKT']
    
    ####################################################################################################################
    etree.indent(root)
    tree = etree.ElementTree(root)
    tree.write(outname, pretty_print=True, xml_declaration=True, encoding='utf-8')


def source_xml(meta, target):
    """
    Function to generate source-level metadata for an NRB target product in OGC 10-157r4 compliant XML format.
    
    Parameters
    ----------
    meta: dict
        Metadata dictionary generated with `metadata.extract.meta_dict`
    target: str
        A path pointing to the root directory of a product scene.
    
    Returns
    -------
    None
    """
    metadir = os.path.join(target, 'source')
    os.makedirs(metadir, exist_ok=True)
    
    for uid in list(meta['source'].keys()):
        scene = os.path.basename(meta['source'][uid]['filename']).split('.')[0]
        outname = os.path.join(metadir, '{}.xml'.format(scene))
        print(outname)
        timeStart = datetime.strftime(meta['source'][uid]['timeStart'], '%Y-%m-%dT%H:%M:%S.%f')
        timeStop = datetime.strftime(meta['source'][uid]['timeStop'], '%Y-%m-%dT%H:%M:%S.%f')
        
        root = etree.Element(_nsc('nrb:EarthObservation'), nsmap=NS_MAP,
                             attrib={_nsc('gml:id'): scene + '_1'})
        
        ################################################################################################################
        phenomenonTime = etree.SubElement(root, _nsc('om:phenomenonTime'))
        timePeriod = etree.SubElement(phenomenonTime, _nsc('gml:TimePeriod'), attrib={_nsc('gml:id'): scene + '_2'})
        beginPosition = etree.SubElement(timePeriod, _nsc('gml:beginPosition'))
        beginPosition.text = timeStart
        endPosition = etree.SubElement(timePeriod, _nsc('gml:endPosition'))
        endPosition.text = timeStop
        
        resultTime = etree.SubElement(root, _nsc('om:resultTime'))
        timeInstant = etree.SubElement(resultTime, _nsc('gml:TimeInstant'), attrib={_nsc('gml:id'): scene + '_3'})
        timePosition = etree.SubElement(timeInstant, _nsc('gml:timePosition'))
        timePosition.text = timeStop
        
        ################################################################################################################
        procedure = etree.SubElement(root, _nsc('om:procedure'))
        earthObservationEquipment = etree.SubElement(procedure, _nsc('eop:EarthObservationEquipment'),
                                                     attrib={_nsc('gml:id'): scene + '_4'})
        platform1, sensor1, acquisition = _common_procedure_elements(eo_equipment=earthObservationEquipment, meta=meta,
                                                                     prod=False)
        
        satReference = etree.SubElement(platform1, _nsc('nrb:satelliteReference'),
                                        attrib={_nsc('xlink:href'): meta['common']['platformReference']})
        radarCenterFreq = etree.SubElement(sensor1, _nsc('nrb:radarCenterFrequency'),
                                           attrib={'uom': 'Hz'})
        radarCenterFreq.text = '{:.3e}'.format(meta['common']['radarCenterFreq'])
        sensorCalibration = etree.SubElement(sensor1, _nsc('nrb:sensorCalibration'),
                                             attrib={_nsc('xlink:href'): meta['source'][uid]['sensorCalibration']})
        
        antennaLookDirection = etree.SubElement(acquisition, _nsc('sar:antennaLookDirection'))
        antennaLookDirection.text = meta['common']['antennaLookDirection']
        orbitMeanAltitude = etree.SubElement(acquisition, _nsc('nrb:orbitMeanAltitude'),
                                             attrib={'uom': 'm'})
        orbitMeanAltitude.text = meta['common']['orbitMeanAltitude']
        orbitDataSource = etree.SubElement(acquisition, _nsc('nrb:orbitDataSource'))
        orbitDataSource.text = meta['source'][uid]['orbitDataSource'].upper()
        ascendingNodeDate = etree.SubElement(acquisition, _nsc('eop:ascendingNodeDate'))
        ascendingNodeDate.text = meta['source'][uid]['ascendingNodeDate']
        startTimeFromAscendingNode = etree.SubElement(acquisition, _nsc('eop:startTimeFromAscendingNode'),
                                                      attrib={'uom': 'ms'})
        startTimeFromAscendingNode.text = meta['source'][uid]['timeStartFromAscendingNode']
        completionTimeFromAscendingNode = etree.SubElement(acquisition, _nsc('eop:completionTimeFromAscendingNode'),
                                                           attrib={'uom': 'ms'})
        completionTimeFromAscendingNode.text = meta['source'][uid]['timeCompletionFromAscendingNode']
        dataTakeID = etree.SubElement(acquisition, _nsc('nrb:dataTakeID'))
        dataTakeID.text = meta['source'][uid]['datatakeID']
        majorCycleID = etree.SubElement(acquisition, _nsc('nrb:majorCycleID'))
        majorCycleID.text = meta['source'][uid]['majorCycleID']
        instrumentAzimuthAngle = etree.SubElement(acquisition, _nsc('eop:instrumentAzimuthAngle'),
                                                  attrib={'uom': 'deg'})
        instrumentAzimuthAngle.text = meta['source'][uid]['instrumentAzimuthAngle']
        minimumIncidenceAngle = etree.SubElement(acquisition, _nsc('sar:minimumIncidenceAngle'),
                                                 attrib={'uom': 'deg'})
        minimumIncidenceAngle.text = str(meta['source'][uid]['incidenceAngleMin'])
        maximumIncidenceAngle = etree.SubElement(acquisition, _nsc('sar:maximumIncidenceAngle'),
                                                 attrib={'uom': 'deg'})
        maximumIncidenceAngle.text = str(meta['source'][uid]['incidenceAngleMax'])
        
        ################################################################################################################
        observedProperty = etree.SubElement(root, _nsc('om:observedProperty'), attrib={'nilReason': 'inapplicable'})
        
        ################################################################################################################
        featureOfInterest = etree.SubElement(root, _nsc('om:featureOfInterest'))
        footprint = etree.SubElement(featureOfInterest, _nsc('eop:Footprint'), attrib={_nsc('gml:id'): scene + '_5'})
        
        multiExtentOf = etree.SubElement(footprint, _nsc('eop:multiExtentOf'))
        multiSurface = etree.SubElement(multiExtentOf, _nsc('gml:MultiSurface'), attrib={_nsc('gml:id'): scene + '_6'})
        surfaceMember = etree.SubElement(multiSurface, _nsc('gml:surfaceMember'))
        polygon = etree.SubElement(surfaceMember, _nsc('gml:Polygon'), attrib={_nsc('gml:id'): scene + '_7'})
        exterior = etree.SubElement(polygon, _nsc('gml:exterior'))
        linearRing = etree.SubElement(exterior, _nsc('gml:LinearRing'))
        posList = etree.SubElement(linearRing, _nsc('gml:posList'), attrib={'uom': 'deg'})
        posList.text = meta['source'][uid]['geom_xml_envelop']
        
        centerOf = etree.SubElement(footprint, _nsc('eop:centerOf'))
        point = etree.SubElement(centerOf, _nsc('gml:Point'), attrib={_nsc('gml:id'): scene + '_8'})
        pos = etree.SubElement(point, _nsc('gml:pos'), attrib={'uom': 'deg'})
        pos.text = meta['source'][uid]['geom_xml_center']
        
        ################################################################################################################
        result = etree.SubElement(root, _nsc('om:result'))
        earthObservationResult = etree.SubElement(result, _nsc('eop:EarthObservationResult'),
                                                  attrib={_nsc('gml:id'): scene + '_9'})
        product = etree.SubElement(earthObservationResult, _nsc('eop:product'))
        productInformation = etree.SubElement(product, _nsc('nrb:ProductInformation'))
        fileName = etree.SubElement(productInformation, _nsc('eop:fileName'))
        serviceReference = etree.SubElement(fileName, _nsc('ows:ServiceReference'), attrib={_nsc('xlink:href'): scene})
        requestMessage = etree.SubElement(serviceReference, _nsc('ows:RequestMessage'))
        
        ################################################################################################################
        metaDataProperty = etree.SubElement(root, _nsc('eop:metaDataProperty'))
        earthObservationMetaData = etree.SubElement(metaDataProperty, _nsc('nrb:EarthObservationMetaData'))
        
        identifier = etree.SubElement(earthObservationMetaData, _nsc('eop:identifier'))
        identifier.text = scene
        doi = etree.SubElement(earthObservationMetaData, _nsc('eop:doi'))
        doi.text = meta['source'][uid]['doi']
        status = etree.SubElement(earthObservationMetaData, _nsc('eop:status'))
        status.text = meta['source'][uid]['status']
        acquisitionType = etree.SubElement(earthObservationMetaData, _nsc('eop:acquisitionType'))
        acquisitionType.text = meta['source'][uid]['acquisitionType']
        productType = etree.SubElement(earthObservationMetaData, _nsc('nrb:productType'),
                                       attrib={'codeSpace': 'urn:esa:eop:Sentinel1:class'})
        productType.text = meta['source'][uid]['productType']
        
        processing = etree.SubElement(earthObservationMetaData, _nsc('eop:processing'))
        processingInformation = etree.SubElement(processing, _nsc('nrb:ProcessingInformation'))
        processingCenter = etree.SubElement(processingInformation, _nsc('eop:processingCenter'),
                                            attrib={'codeSpace': 'urn:esa:eop:Sentinel1:facility'})
        processingCenter.text = meta['source'][uid]['processingCenter']
        processingDate = etree.SubElement(processingInformation, _nsc('eop:processingDate'))
        processingDate.text = meta['source'][uid]['processingDate']
        processorName = etree.SubElement(processingInformation, _nsc('eop:processorName'))
        processorName.text = meta['source'][uid]['processorName']
        processorVersion = etree.SubElement(processingInformation, _nsc('eop:processorVersion'))
        processorVersion.text = meta['source'][uid]['processorVersion']
        processingLevel = etree.SubElement(processingInformation, _nsc('eop:processingLevel'))
        processingLevel.text = meta['common']['processingLevel']
        processingMode = etree.SubElement(processingInformation, _nsc('eop:processingMode'))
        processingMode.text = meta['source'][uid]['processingMode']
        orbitStateVector = etree.SubElement(processingInformation, _nsc('nrb:orbitStateVector'),
                                            attrib={'access': meta['source'][uid]['orbitDataAccess']})
        orbitStateVector.text = meta['source'][uid]['orbitStateVector']
        for swath in meta['source'][uid]['swaths']:
            azimuthLookBandwidth = etree.SubElement(processingInformation, _nsc('nrb:azimuthLookBandwidth'),
                                                    attrib={'uom': 'Hz', 'beam': swath})
            azimuthLookBandwidth.text = str(meta['source'][uid]['azimuthLookBandwidth'][swath])
        for swath in meta['source'][uid]['swaths']:
            rangeLookBandwidth = etree.SubElement(processingInformation, _nsc('nrb:rangeLookBandwidth'),
                                                  attrib={'uom': 'Hz', 'beam': swath})
            rangeLookBandwidth.text = str(meta['source'][uid]['rangeLookBandwidth'][swath])
        lutApplied = etree.SubElement(processingInformation, _nsc('nrb:lutApplied'))
        lutApplied.text = meta['source'][uid]['lutApplied']
        
        performance = etree.SubElement(earthObservationMetaData, _nsc('nrb:performance'))
        performanceIndicators = etree.SubElement(performance, _nsc('nrb:PerformanceIndicators'))
        noiseEquivalentIntensityType = etree.SubElement(performanceIndicators, _nsc('nrb:noiseEquivalentIntensityType'),
                                                        attrib={'uom': 'dB'})
        noiseEquivalentIntensityType.text = str(meta['source'][uid]['perfNoiseEquivalentIntensityType'])
        for pol in meta['common']['polarisationChannels']:
            estimatesMin = etree.SubElement(performanceIndicators, _nsc('nrb:estimates'),
                                            attrib={'pol': pol, 'type': 'minimum'})
            estimatesMin.text = str(meta['source'][uid]['perfEstimates'][pol]['minimum'])
            estimatesMax = etree.SubElement(performanceIndicators, _nsc('nrb:estimates'),
                                            attrib={'pol': pol, 'type': 'maximum'})
            estimatesMax.text = str(meta['source'][uid]['perfEstimates'][pol]['maximum'])
            estimatesMean = etree.SubElement(performanceIndicators, _nsc('nrb:estimates'),
                                             attrib={'pol': pol, 'type': 'mean'})
            estimatesMean.text = str(meta['source'][uid]['perfEstimates'][pol]['mean'])
        equivalentNumberOfLooks = etree.SubElement(performanceIndicators, _nsc('nrb:equivalentNumberOfLooks'))
        equivalentNumberOfLooks.text = str(meta['source'][uid]['perfEquivalentNumberOfLooks'])
        peakSideLobeRatio = etree.SubElement(performanceIndicators, _nsc('nrb:peakSideLobeRatio'),
                                             attrib={'uom': 'dB'})
        peakSideLobeRatio.text = str(meta['source'][uid]['perfPeakSideLobeRatio'])
        integratedSideLobeRatio = etree.SubElement(performanceIndicators, _nsc('nrb:integratedSideLobeRatio'),
                                                   attrib={'uom': 'dB'})
        integratedSideLobeRatio.text = str(meta['source'][uid]['perfIntegratedSideLobeRatio'])
        
        azimuthNumberOfLooks = etree.SubElement(earthObservationMetaData, _nsc('nrb:azimuthNumberOfLooks'))
        azimuthNumberOfLooks.text = meta['source'][uid]['azimuthNumberOfLooks']
        rangeNumberOfLooks = etree.SubElement(earthObservationMetaData, _nsc('nrb:rangeNumberOfLooks'))
        rangeNumberOfLooks.text = meta['source'][uid]['rangeNumberOfLooks']
        dataGeometry = etree.SubElement(earthObservationMetaData, _nsc('nrb:dataGeometry'))
        dataGeometry.text = meta['source'][uid]['dataGeometry']
        for swath in meta['source'][uid]['swaths']:
            azimuthResolution = etree.SubElement(earthObservationMetaData, _nsc('nrb:azimuthResolution'),
                                                 attrib={'uom': 'm', 'beam': swath})
            azimuthResolution.text = meta['source'][uid]['azimuthResolution'][swath]
        for swath in meta['source'][uid]['swaths']:
            rangeResolution = etree.SubElement(earthObservationMetaData, _nsc('nrb:rangeResolution'),
                                               attrib={'uom': 'm', 'beam': swath})
            rangeResolution.text = meta['source'][uid]['rangeResolution'][swath]
        azimuthPixelSpacing = etree.SubElement(earthObservationMetaData, _nsc('nrb:azimuthPixelSpacing'),
                                               attrib={'uom': 'm'})
        azimuthPixelSpacing.text = meta['source'][uid]['azimuthPixelSpacing']
        rangePixelSpacing = etree.SubElement(earthObservationMetaData, _nsc('nrb:rangePixelSpacing'),
                                             attrib={'uom': 'm'})
        rangePixelSpacing.text = meta['source'][uid]['rangePixelSpacing']
        polCalMatrices = etree.SubElement(earthObservationMetaData, _nsc('nrb:polCalMatrices'),
                                          attrib={_nsc('xlink:href'): str(meta['source'][uid]['polCalMatrices'])})
        meanFaradayRotationAngle = etree.SubElement(earthObservationMetaData, _nsc('nrb:meanFaradayRotationAngle'),
                                                    attrib={'uom': 'deg'})
        meanFaradayRotationAngle.text = meta['source'][uid]['faradayMeanRotationAngle']
        faraday_ref = str(meta['source'][uid]['faradayRotationReference'])
        referenceFaradayRotation = etree.SubElement(earthObservationMetaData, _nsc('nrb:referenceFaradayRotation'),
                                                    attrib={_nsc('xlink:href'): faraday_ref})
        ionosphereIndicator = etree.SubElement(earthObservationMetaData, _nsc('nrb:ionosphereIndicator'))
        ionosphereIndicator.text = meta['source'][uid]['ionosphereIndicator']
        
        ################################################################################################################
        etree.indent(root)
        tree = etree.ElementTree(root)
        tree.write(outname, pretty_print=True, xml_declaration=True, encoding='utf-8')


def main(meta, target, tifs):
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
    
    Returns
    -------
    None
    """
    source_xml(meta=meta, target=target)
    product_xml(meta=meta, target=target, tifs=tifs)
