import os
import re
from lxml import etree
from datetime import datetime
from spatialist import Raster
from S1_NRB.metadata.mapping import SAMPLE_MAP


nsmap_out = {'nrb': 'http://earth.esa.int/sentinel-1/nrb/1.0',
             'eop': 'http://www.opengis.net/eop/2.1',
             'gml': 'http://www.opengis.net/gml/3.2',
             'om': 'http://www.opengis.net/om/2.0',
             'ows': 'http://www.opengis.net/ows/2.0',
             'sar': 'http://www.opengis.net/sar/2.1',
             'xlink': 'http://www.w3.org/1999/xlink',
             'xsi': 'http://www.w3.org/2001/XMLSchema-instance'}


def _nsc(text):
    ns, key = text.split(':')
    return '{{{0}}}{1}'.format(nsmap_out[ns], key)


def _get_ref_type(ref_link):
    if ref_link is None:
        return 'None'
    elif 'doi.org' in ref_link:
        return 'DOI'
    else:
        return 'URL'


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
    
    root = etree.Element(_nsc('nrb:EarthObservation'), nsmap=nsmap_out,
                         attrib={_nsc('gml:id'): scene_id + '_1'})
    
    ####################################################################################################################
    phenomenonTime = etree.SubElement(root, _nsc('om:phenomenonTime'))
    TimePeriod = etree.SubElement(phenomenonTime, _nsc('gml:TimePeriod'),
                                  attrib={_nsc('gml:id'): scene_id + '_2'})
    beginPosition = etree.SubElement(TimePeriod, _nsc('gml:beginPosition'))
    beginPosition.text = timeStart
    endPosition = etree.SubElement(TimePeriod, _nsc('gml:endPosition'))
    endPosition.text = timeStop
    
    ####################################################################################################################
    resultTime = etree.SubElement(root, _nsc('om:resultTime'))
    TimeInstant = etree.SubElement(resultTime, _nsc('gml:TimeInstant'),
                                   attrib={_nsc('gml:id'): scene_id + '_3'})
    timePosition = etree.SubElement(TimeInstant, _nsc('gml:timePosition'))
    timePosition.text = timeStop
    
    ####################################################################################################################
    procedure = etree.SubElement(root, _nsc('om:procedure'))
    EarthObservationEquipment = etree.SubElement(procedure, _nsc('eop:EarthObservationEquipment'),
                                                 attrib={_nsc('gml:id'): scene_id + '_4'})
    
    platform = etree.SubElement(EarthObservationEquipment, _nsc('eop:platform'))
    Platform = etree.SubElement(platform, _nsc('eop:Platform'))
    shortName = etree.SubElement(Platform, _nsc('eop:shortName'))
    shortName.text = meta['common']['platformShortName'].upper()
    serialIdentifier = etree.SubElement(Platform, _nsc('eop:serialIdentifier'))
    serialIdentifier.text = meta['common']['platformIdentifier']
    
    instrument = etree.SubElement(EarthObservationEquipment, _nsc('eop:instrument'))
    Instrument = etree.SubElement(instrument, _nsc('eop:Instrument'))
    shortName = etree.SubElement(Instrument, _nsc('eop:shortName'))
    shortName.text = meta['common']['instrumentShortName']
    
    sensor = etree.SubElement(EarthObservationEquipment, _nsc('eop:sensor'))
    Sensor = etree.SubElement(sensor, _nsc('nrb:Sensor'))
    sensorType = etree.SubElement(Sensor, _nsc('eop:sensorType'))
    sensorType.text = meta['common']['sensorType']
    radarBand = etree.SubElement(Sensor, _nsc('nrb:radarBand'))
    radarBand.text = meta['common']['radarBand']
    operationalMode = etree.SubElement(Sensor, _nsc('eop:operationalMode'),
                                       attrib={'codeSpace': 'urn:esa:eop:C-SAR:operationalMode'})
    operationalMode.text = meta['common']['operationalMode']
    swathIdentifier = etree.SubElement(Sensor, _nsc('eop:swathIdentifier'),
                                       attrib={'codeSpace': 'urn:esa:eop:C-SAR:swathIdentifier'})
    swathIdentifier.text = meta['common']['swathIdentifier']
    
    acquisitionParameters = etree.SubElement(EarthObservationEquipment, _nsc('eop:acquisitionParameters'))
    Acquisition = etree.SubElement(acquisitionParameters, _nsc('nrb:Acquisition'))
    polarisationMode = etree.SubElement(Acquisition, _nsc('sar:polarisationMode'))
    polarisationMode.text = meta['common']['polarisationMode']
    polarisationChannels = etree.SubElement(Acquisition, _nsc('sar:polarisationChannels'))
    polarisationChannels.text = ', '.join(meta['common']['polarisationChannels'])
    orbitDirection = etree.SubElement(Acquisition, _nsc('eop:orbitDirection'))
    orbitDirection.text = meta['common']['orbitDirection'].upper()
    orbitNumber = etree.SubElement(Acquisition, _nsc('eop:orbitNumber'))
    orbitNumber.text = meta['common']['orbitNumber']
    wrsLongitudeGrid = etree.SubElement(Acquisition, _nsc('eop:wrsLongitudeGrid'),
                                        attrib={'codeSpace': 'urn:esa:eop:Sentinel1:relativeOrbits'})
    wrsLongitudeGrid.text = meta['common']['wrsLongitudeGrid']
    numberOfAcquisitions = etree.SubElement(Acquisition, _nsc('nrb:numberOfAcquisitions'))
    numberOfAcquisitions.text = meta['prod']['numberOfAcquisitions']
    
    ####################################################################################################################
    observedProperty = etree.SubElement(root, _nsc('om:observedProperty'),
                                        attrib={_nsc('xsi:nil'): 'true', 'nilReason': 'inapplicable'})
    
    ####################################################################################################################
    featureOfInterest = etree.SubElement(root, _nsc('om:featureOfInterest'))
    Footprint = etree.SubElement(featureOfInterest, _nsc('eop:Footprint'), attrib={_nsc('gml:id'): scene_id + '_5'})
    
    multiExtentOf = etree.SubElement(Footprint, _nsc('eop:multiExtentOf'))
    MultiSurface = etree.SubElement(multiExtentOf, _nsc('gml:MultiSurface'), attrib={_nsc('gml:id'): scene_id + '_6'})
    surfaceMember = etree.SubElement(MultiSurface, _nsc('gml:surfaceMember'))
    Polygon = etree.SubElement(surfaceMember, _nsc('gml:Polygon'), attrib={_nsc('gml:id'): scene_id + '_7'})
    exterior = etree.SubElement(Polygon, _nsc('gml:exterior'))
    LinearRing = etree.SubElement(exterior, _nsc('gml:LinearRing'))
    posList = etree.SubElement(LinearRing, _nsc('gml:posList'), attrib={'uom': 'deg'})
    posList.text = meta['prod']['geom_xml_envelope']
    
    centerOf = etree.SubElement(Footprint, _nsc('eop:centerOf'))
    Point = etree.SubElement(centerOf, _nsc('gml:Point'), attrib={_nsc('gml:id'): scene_id + '_8'})
    pos = etree.SubElement(Point, _nsc('gml:pos'), attrib={'uom': 'deg'})
    pos.text = meta['prod']['geom_xml_center']
    
    ####################################################################################################################
    result = etree.SubElement(root, _nsc('om:result'))
    EarthObservationResult = etree.SubElement(result, _nsc('eop:EarthObservationResult'),
                                              attrib={_nsc('gml:id'): scene_id + '_9'})
    
    product = etree.SubElement(EarthObservationResult, _nsc('eop:product'))
    ProductInformation = etree.SubElement(product, _nsc('nrb:ProductInformation'))
    fileName = etree.SubElement(ProductInformation, _nsc('eop:fileName'))
    ServiceReference = etree.SubElement(fileName, _nsc('ows:ServiceReference'),
                                        attrib={_nsc('xlink:href'): scene_id})
    RequestMessage = etree.SubElement(ServiceReference, _nsc('ows:RequestMessage'))
    version = etree.SubElement(ProductInformation, _nsc('eop:size'))
    
    for tif in tifs:
        relpath = './' + os.path.relpath(tif, target).replace('\\', '/')
        product = etree.SubElement(EarthObservationResult, _nsc('eop:product'))
        ProductInformation = etree.SubElement(product, _nsc('nrb:ProductInformation'))
        fileName = etree.SubElement(ProductInformation, _nsc('eop:fileName'))
        ServiceReference = etree.SubElement(fileName, _nsc('ows:ServiceReference'),
                                            attrib={_nsc('xlink:href'): relpath})
        RequestMessage = etree.SubElement(ServiceReference, _nsc('ows:RequestMessage'))
        size = etree.SubElement(ProductInformation, _nsc('eop:size'), attrib={'uom': 'bytes'})
        size.text = str(os.path.getsize(tif))
        byteOrder = etree.SubElement(ProductInformation, _nsc('nrb:byteOrder'))
        byteOrder.text = meta['prod']['fileByteOrder']
        dataFormat = etree.SubElement(ProductInformation, _nsc('nrb:dataFormat'))
        dataFormat.text = meta['prod']['fileFormat']
        dataType = etree.SubElement(ProductInformation, _nsc('nrb:dataType'))
        dataType.text = meta['prod']['fileDataType'].upper()
        bitsPerSample = etree.SubElement(ProductInformation, _nsc('nrb:bitsPerSample'))
        bitsPerSample.text = meta['prod']['fileBitsPerSample']
        noDataVal = etree.SubElement(ProductInformation, _nsc('nrb:noDataValue'))
        noDataVal.text = 'NaN'
        
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
                        bands = dm_ras.bands
                    if bands > 1:   # multi-band data mask (default)
                        samples = list(SAMPLE_MAP[key]['values'].values())
                        samples.remove('layover and shadow')
                        if bands != len(samples):
                            raise RuntimeError('Mismatch between number of bands ({nbands}) of the '
                                               'multi-band data mask file and the number of keys '
                                               'in SAMPLE_MAP ({nkeys}).'.format(nbands=bands, nkeys=len(samples)))
                        for i in range(bands):
                            bitValue = etree.SubElement(ProductInformation, _nsc('nrb:bitValue'),
                                                        attrib={'band': str(i+1), 'name': samples[i]})
                            bitValue.text = '1'
                        pass
                    else:  # single-band data mask
                        for val in SAMPLE_MAP[key]['values']:
                            bitValue = etree.SubElement(ProductInformation, _nsc('nrb:bitValue'),
                                                        attrib={'band': '1', 'name': SAMPLE_MAP[key]['values'][val]})
                            bitValue.text = str(val)
                else:  # key == '-id.tif'
                    src_list = list(meta['source'].keys())
                    src_target = [os.path.basename(meta['source'][src]['filename']).replace('.SAFE', '').replace('.zip', '')
                                  for src in src_list]
                    for i, s in enumerate(src_target):
                        bitValue = etree.SubElement(ProductInformation, _nsc('nrb:bitValue'),
                                                    attrib={'band': '1', 'name': s})
                        bitValue.text = str(i+1)
            
            if key == '-ei.tif':
                ellipsoidalHeight = etree.SubElement(ProductInformation, _nsc('nrb:ellipsoidalHeight'),
                                                     attrib={'uom': 'm'})
                ellipsoidalHeight.text = meta['prod']['ellipsoidalHeight']
            
            if SAMPLE_MAP[key]['unit'] is None:
                sampleType = etree.SubElement(ProductInformation, _nsc('nrb:sampleType'))
            else:
                sampleType = etree.SubElement(ProductInformation, _nsc('nrb:sampleType'),
                                              attrib={'uom': SAMPLE_MAP[key]['unit']})
            sampleType.text = SAMPLE_MAP[key]['type']
        
        if 'measurement' in tif:
            creationTime = etree.SubElement(ProductInformation, _nsc('nrb:creationTime'))
            creationTime.text = datetime.fromtimestamp(os.path.getctime(tif)).isoformat()
            polarization = etree.SubElement(ProductInformation, _nsc('nrb:polarization'))
            polarization.text = re.search('[vh]{2}', tif).group().upper()
            numBorderPixels = etree.SubElement(ProductInformation, _nsc('nrb:numBorderPixels'))
            numBorderPixels.text = str(meta['prod']['numBorderPixels'])
    
    ####################################################################################################################
    metaDataProperty = etree.SubElement(root, _nsc('eop:metaDataProperty'))
    EarthObservationMetaData = etree.SubElement(metaDataProperty, _nsc('nrb:EarthObservationMetaData'))
    
    identifier = etree.SubElement(EarthObservationMetaData, _nsc('eop:identifier'))
    identifier.text = scene_id
    doi = etree.SubElement(EarthObservationMetaData, _nsc('eop:doi'))
    doi.text = meta['prod']['doi']
    status = etree.SubElement(EarthObservationMetaData, _nsc('eop:status'))
    status.text = meta['prod']['status']
    acquisitionType = etree.SubElement(EarthObservationMetaData, _nsc('eop:acquisitionType'))
    acquisitionType.text = meta['prod']['acquisitionType']
    productType = etree.SubElement(EarthObservationMetaData, _nsc('nrb:productType'),
                                   attrib={'codeSpace': 'urn:esa:eop:Sentinel1:class'})
    productType.text = meta['prod']['productName-short']
    refDoc = etree.SubElement(EarthObservationMetaData, _nsc('nrb:refDoc'),
                              attrib={'name': meta['prod']['productName'],
                                      'version': meta['prod']['card4l-version'],
                                      _nsc('xlink:href'): meta['prod']['card4l-link']})
    
    processing = etree.SubElement(EarthObservationMetaData, _nsc('eop:processing'))
    ProcessingInformation = etree.SubElement(processing, _nsc('nrb:ProcessingInformation'))
    processingCenter = etree.SubElement(ProcessingInformation, _nsc('eop:processingCenter'),
                                        attrib={'codeSpace': 'urn:esa:eop:Sentinel1:facility'})
    processingCenter.text = meta['prod']['processingCenter']
    processingDate = etree.SubElement(ProcessingInformation, _nsc('eop:processingDate'))
    processingDate.text = timeCreated
    processorName = etree.SubElement(ProcessingInformation, _nsc('eop:processorName'))
    processorName.text = meta['prod']['processorName']
    processorVersion = etree.SubElement(ProcessingInformation, _nsc('eop:processorVersion'))
    processorVersion.text = meta['prod']['processorVersion']
    processingLevel = etree.SubElement(ProcessingInformation, _nsc('eop:processingLevel'))
    processingLevel.text = meta['common']['processingLevel']
    processingMode = etree.SubElement(ProcessingInformation, _nsc('eop:processingMode'),
                                      attrib={'codeSpace': 'urn:esa:eop:Sentinel1:class'})
    processingMode.text = meta['prod']['processingMode']
    for src in list(meta['source'].keys()):
        src_path = '{}.xml'.format(os.path.basename(meta['source'][src]['filename']).split('.')[0])
        src_target = os.path.join('./source', src_path).replace('\\', '/')
        sourceProduct = etree.SubElement(ProcessingInformation, _nsc('nrb:sourceProduct'),
                                         attrib={_nsc('xlink:href'): src_target})
    speckleFilterApplied = etree.SubElement(ProcessingInformation, _nsc('nrb:speckleFilterApplied'))
    speckleFilterApplied.text = str(meta['prod']['speckleFilterApplied']).lower()
    NRApplied = etree.SubElement(ProcessingInformation, _nsc('nrb:NRApplied'))
    NRApplied.text = str(meta['prod']['NRApplied']).lower()
    if meta['prod']['NRApplied']:
        NRAlgorithm = etree.SubElement(ProcessingInformation, _nsc('nrb:NRAlgorithm'),
                                       attrib={_nsc('xlink:href'): meta['prod']['NRAlgorithm']})
    RTCAlgorithm = etree.SubElement(ProcessingInformation, _nsc('nrb:RTCAlgorithm'),
                                    attrib={_nsc('xlink:href'): meta['prod']['RTCAlgorithm']})
    geoCorrAlgorithm = etree.SubElement(ProcessingInformation, _nsc('nrb:geoCorrAlgorithm'),
                                        attrib={_nsc('xlink:href'): meta['prod']['geoCorrAlgorithm']})
    geoCorrResamplingMethod = etree.SubElement(ProcessingInformation, _nsc('nrb:geoCorrResamplingAlgorithm'))
    geoCorrResamplingMethod.text = meta['prod']['geoCorrResamplingMethod'].upper()
    DEMReference = etree.SubElement(ProcessingInformation, _nsc('nrb:DEMReference'),
                                    attrib={'name': meta['prod']['demName'],
                                            'dem': meta['prod']['demType'],
                                            _nsc('xlink:href'): meta['prod']['demReference']})
    DEMResamplingMethod = etree.SubElement(ProcessingInformation, _nsc('nrb:DEMResamplingMethod'))
    DEMResamplingMethod.text = meta['prod']['demResamplingMethod'].upper()
    DEMAccess = etree.SubElement(ProcessingInformation, _nsc('nrb:DEMAccess'),
                                 attrib={_nsc('xlink:href'): meta['prod']['demAccess']})
    EGMReference = etree.SubElement(ProcessingInformation, _nsc('nrb:EGMReference'),
                                    attrib={_nsc('xlink:href'): meta['prod']['demEGMReference']})
    egmResamplingMethod = etree.SubElement(ProcessingInformation, _nsc('nrb:EGMResamplingMethod'))
    egmResamplingMethod.text = meta['prod']['demEGMResamplingMethod'].upper()
    
    radiometricAccuracyRelative = etree.SubElement(EarthObservationMetaData, _nsc('nrb:radiometricAccuracyRelative'),
                                                   attrib={'uom': 'dB'})
    radiometricAccuracyRelative.text = meta['prod']['radiometricAccuracyRelative']
    radiometricAccuracyAbsolute = etree.SubElement(EarthObservationMetaData, _nsc('nrb:radiometricAccuracyAbsolute'),
                                                   attrib={'uom': 'dB'})
    radiometricAccuracyAbsolute.text = meta['prod']['radiometricAccuracyAbsolute']
    radiometricAccuracyReference = etree.SubElement(EarthObservationMetaData, _nsc('nrb:radiometricAccuracyReference'),
                                                    attrib={_nsc('xlink:href'): str(meta['prod']['radiometricAccuracyReference'])})
    geoCorrAccuracyType = etree.SubElement(EarthObservationMetaData, _nsc('nrb:geoCorrAccuracyType'))
    geoCorrAccuracyType.text = meta['prod']['geoCorrAccuracyType']
    geoCorrAccuracyNorthernSTDev = etree.SubElement(EarthObservationMetaData, _nsc('nrb:geoCorrAccuracyNorthernSTDev'),
                                                    attrib={'uom': 'm'})
    geoCorrAccuracyNorthernSTDev.text = meta['prod']['geoCorrAccuracyNorthernSTDev']
    geoCorrAccuracyEasternSTDev = etree.SubElement(EarthObservationMetaData, _nsc('nrb:geoCorrAccuracyEasternSTDev'),
                                                   attrib={'uom': 'm'})
    geoCorrAccuracyEasternSTDev.text = meta['prod']['geoCorrAccuracyEasternSTDev']
    geoCorrAccuracyNorthernBias = etree.SubElement(EarthObservationMetaData, _nsc('nrb:geoCorrAccuracyNorthernBias'),
                                                   attrib={'uom': 'm'})
    geoCorrAccuracyNorthernBias.text = meta['prod']['geoCorrAccuracyNorthernBias']
    geoCorrAccuracyEasternBias = etree.SubElement(EarthObservationMetaData, _nsc('nrb:geoCorrAccuracyEasternBias'),
                                                  attrib={'uom': 'm'})
    geoCorrAccuracyEasternBias.text = meta['prod']['geoCorrAccuracyEasternBias']
    geoCorrAccuracy_rRMSE = etree.SubElement(EarthObservationMetaData, _nsc('nrb:geoCorrAccuracy_rRMSE'),
                                             attrib={'uom': 'm'})
    geoCorrAccuracy_rRMSE.text = meta['prod']['geoCorrAccuracy_rRMSE']
    geoCorrAccuracyReference = etree.SubElement(EarthObservationMetaData, _nsc('nrb:geoCorrAccuracyReference'),
                                                attrib={_nsc('xlink:href'): meta['prod']['geoCorrAccuracyReference']})
    azimuthNumberOfLooks = etree.SubElement(EarthObservationMetaData, _nsc('nrb:azimuthNumberOfLooks'))
    azimuthNumberOfLooks.text = meta['prod']['azimuthNumberOfLooks']
    rangeNumberOfLooks = etree.SubElement(EarthObservationMetaData, _nsc('nrb:rangeNumberOfLooks'))
    rangeNumberOfLooks.text = meta['prod']['rangeNumberOfLooks']
    numLines = etree.SubElement(EarthObservationMetaData, _nsc('nrb:numLines'))
    numLines.text = meta['prod']['numLines']
    numPixelsPerLine = etree.SubElement(EarthObservationMetaData, _nsc('nrb:numPixelsPerLine'))
    numPixelsPerLine.text = meta['prod']['numPixelsPerLine']
    columnSpacing = etree.SubElement(EarthObservationMetaData, _nsc('nrb:columnSpacing'), attrib={'uom': 'm'})
    columnSpacing.text = meta['prod']['pxSpacingColumn']
    rowSpacing = etree.SubElement(EarthObservationMetaData, _nsc('nrb:rowSpacing'), attrib={'uom': 'm'})
    rowSpacing.text = meta['prod']['pxSpacingRow']
    pixelCoordinateConvention = etree.SubElement(EarthObservationMetaData, _nsc('nrb:pixelCoordinateConvention'))
    pixelCoordinateConvention.text = meta['prod']['pixelCoordinateConvention']
    backscatterMeasurement = etree.SubElement(EarthObservationMetaData, _nsc('nrb:backscatterMeasurement'))
    backscatterMeasurement.text = meta['prod']['backscatterMeasurement']
    backscatterConvention = etree.SubElement(EarthObservationMetaData, _nsc('nrb:backscatterConvention'))
    backscatterConvention.text = meta['prod']['backscatterConvention']
    backscatterConversionEq = etree.SubElement(EarthObservationMetaData, _nsc('nrb:backscatterConversionEq'),
                                               attrib={'uom': 'dB'})
    backscatterConversionEq.text = meta['prod']['backscatterConversionEq']
    griddingConvention = etree.SubElement(EarthObservationMetaData, _nsc('nrb:griddingConvention'),
                                          attrib={_nsc('xlink:href'): meta['prod']['griddingConventionURL']})
    mgrsID = etree.SubElement(EarthObservationMetaData, _nsc('nrb:mgrsID'))
    mgrsID.text = meta['prod']['mgrsID']
    crsEPSG = etree.SubElement(EarthObservationMetaData, _nsc('nrb:crsEPSG'),
                               attrib={'codespace': 'urn:esa:eop:crs'})
    crsEPSG.text = meta['prod']['crsEPSG']
    crsWKT = etree.SubElement(EarthObservationMetaData, _nsc('nrb:crsWKT'))
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
        
        root = etree.Element(_nsc('nrb:EarthObservation'), nsmap=nsmap_out,
                             attrib={_nsc('gml:id'): scene + '_1'})
        
        ################################################################################################################
        phenomenonTime = etree.SubElement(root, _nsc('om:phenomenonTime'))
        TimePeriod = etree.SubElement(phenomenonTime, _nsc('gml:TimePeriod'),
                                      attrib={_nsc('gml:id'): scene + '_2'})
        beginPosition = etree.SubElement(TimePeriod, _nsc('gml:beginPosition'))
        beginPosition.text = timeStart
        endPosition = etree.SubElement(TimePeriod, _nsc('gml:endPosition'))
        endPosition.text = timeStop
        
        ################################################################################################################
        resultTime = etree.SubElement(root, _nsc('om:resultTime'))
        TimeInstant = etree.SubElement(resultTime, _nsc('gml:TimeInstant'),
                                       attrib={_nsc('gml:id'): scene + '_3'})
        timePosition = etree.SubElement(TimeInstant, _nsc('gml:timePosition'))
        timePosition.text = timeStop
        
        ################################################################################################################
        procedure = etree.SubElement(root, _nsc('om:procedure'))
        EarthObservationEquipment = etree.SubElement(procedure, _nsc('eop:EarthObservationEquipment'),
                                                     attrib={_nsc('gml:id'): scene + '_4'})
        
        platform = etree.SubElement(EarthObservationEquipment, _nsc('eop:platform'))
        Platform = etree.SubElement(platform, _nsc('eop:Platform'))
        shortName = etree.SubElement(Platform, _nsc('eop:shortName'))
        shortName.text = meta['common']['platformShortName'].upper()
        serialIdentifier = etree.SubElement(Platform, _nsc('eop:serialIdentifier'))
        serialIdentifier.text = meta['common']['platformIdentifier']
        satReference = etree.SubElement(Platform, _nsc('nrb:satelliteReference'),
                                        attrib={_nsc('xlink:href'): meta['common']['platformReference']})
        
        instrument = etree.SubElement(EarthObservationEquipment, _nsc('eop:instrument'))
        Instrument = etree.SubElement(instrument, _nsc('eop:Instrument'))
        shortName = etree.SubElement(Instrument, _nsc('eop:shortName'))
        shortName.text = meta['common']['instrumentShortName']
        
        sensor = etree.SubElement(EarthObservationEquipment, _nsc('eop:sensor'))
        Sensor = etree.SubElement(sensor, _nsc('nrb:Sensor'))
        sensorType = etree.SubElement(Sensor, _nsc('eop:sensorType'))
        sensorType.text = meta['common']['sensorType']
        radarBand = etree.SubElement(Sensor, _nsc('nrb:radarBand'))
        radarBand.text = meta['common']['radarBand']
        operationalMode = etree.SubElement(Sensor, _nsc('eop:operationalMode'),
                                           attrib={'codeSpace': 'urn:esa:eop:C-SAR:operationalMode'})
        operationalMode.text = meta['common']['operationalMode']
        swathIdentifier = etree.SubElement(Sensor, _nsc('eop:swathIdentifier'),
                                           attrib={'codeSpace': 'urn:esa:eop:C-SAR:swathIdentifier'})
        swathIdentifier.text = meta['common']['swathIdentifier']
        radarCenterFreq = etree.SubElement(Sensor, _nsc('nrb:radarCenterFrequency'),
                                           attrib={'uom': 'Hz'})
        radarCenterFreq.text = '{:.3e}'.format(meta['common']['radarCenterFreq'])
        sensorCalibration = etree.SubElement(Sensor, _nsc('nrb:sensorCalibration'),
                                             attrib={_nsc('xlink:href'): meta['source'][uid]['sensorCalibration']})
        sensorCalibration.text = meta['source'][uid]['sensorCalibration']
        
        acquisitionParameters = etree.SubElement(EarthObservationEquipment, _nsc('eop:acquisitionParameters'))
        Acquisition = etree.SubElement(acquisitionParameters, _nsc('nrb:Acquisition'))
        polarisationMode = etree.SubElement(Acquisition, _nsc('sar:polarisationMode'))
        polarisationMode.text = meta['common']['polarisationMode']
        polarisationChannels = etree.SubElement(Acquisition, _nsc('sar:polarisationChannels'))
        polarisationChannels.text = ', '.join(meta['common']['polarisationChannels'])
        orbitDirection = etree.SubElement(Acquisition, _nsc('eop:orbitDirection'))
        orbitDirection.text = meta['common']['orbitDirection'].upper()
        orbitNumber = etree.SubElement(Acquisition, _nsc('eop:orbitNumber'))
        orbitNumber.text = meta['common']['orbitNumber']
        wrsLongitudeGrid = etree.SubElement(Acquisition, _nsc('eop:wrsLongitudeGrid'),
                                            attrib={'codeSpace': 'urn:esa:eop:Sentinel1:relativeOrbits'})
        wrsLongitudeGrid.text = meta['common']['wrsLongitudeGrid']
        antennaLookDirection = etree.SubElement(Acquisition, _nsc('sar:antennaLookDirection'))
        antennaLookDirection.text = meta['common']['antennaLookDirection']
        orbitMeanAltitude = etree.SubElement(Acquisition, _nsc('nrb:orbitMeanAltitude'),
                                             attrib={'uom': 'm'})
        orbitMeanAltitude.text = meta['common']['orbitMeanAltitude']
        orbitDataSource = etree.SubElement(Acquisition, _nsc('eop:orbitDataSource'))
        orbitDataSource.text = meta['source'][uid]['orbitDataSource'].upper()
        ascendingNodeDate = etree.SubElement(Acquisition, _nsc('eop:ascendingNodeDate'))
        ascendingNodeDate.text = meta['source'][uid]['ascendingNodeDate']
        startTimeFromAscendingNode = etree.SubElement(Acquisition, _nsc('eop:startTimeFromAscendingNode'),
                                                      attrib={'uom': 'ms'})
        startTimeFromAscendingNode.text = meta['source'][uid]['timeStartFromAscendingNode']
        completionTimeFromAscendingNode = etree.SubElement(Acquisition, _nsc('eop:completionTimeFromAscendingNode'),
                                                           attrib={'uom': 'ms'})
        completionTimeFromAscendingNode.text = meta['source'][uid]['timeCompletionFromAscendingNode']
        dataTakeID = etree.SubElement(Acquisition, _nsc('nrb:dataTakeID'))
        dataTakeID.text = meta['source'][uid]['datatakeID']
        majorCycleID = etree.SubElement(Acquisition, _nsc('nrb:majorCycleID'))
        majorCycleID.text = meta['source'][uid]['majorCycleID']
        instrumentAzimuthAngle = etree.SubElement(Acquisition, _nsc('eop:instrumentAzimuthAngle'),
                                                  attrib={'uom': 'deg'})
        instrumentAzimuthAngle.text = meta['source'][uid]['instrumentAzimuthAngle']
        minimumIncidenceAngle = etree.SubElement(Acquisition, _nsc('sar:minimumIncidenceAngle'),
                                                 attrib={'uom': 'deg'})
        minimumIncidenceAngle.text = str(meta['source'][uid]['incidenceAngleMin'])
        maximumIncidenceAngle = etree.SubElement(Acquisition, _nsc('sar:maximumIncidenceAngle'),
                                                 attrib={'uom': 'deg'})
        maximumIncidenceAngle.text = str(meta['source'][uid]['incidenceAngleMax'])
        
        ################################################################################################################
        observedProperty = etree.SubElement(root, _nsc('om:observedProperty'),
                                            attrib={_nsc('xsi:nil'): 'true', 'nilReason': 'inapplicable'})
        
        ################################################################################################################
        featureOfInterest = etree.SubElement(root, _nsc('om:featureOfInterest'))
        Footprint = etree.SubElement(featureOfInterest, _nsc('eop:Footprint'), attrib={_nsc('gml:id'): scene + '_5'})
        
        multiExtentOf = etree.SubElement(Footprint, _nsc('eop:multiExtentOf'))
        MultiSurface = etree.SubElement(multiExtentOf, _nsc('gml:MultiSurface'), attrib={_nsc('gml:id'): scene + '_6'})
        surfaceMember = etree.SubElement(MultiSurface, _nsc('gml:surfaceMember'))
        Polygon = etree.SubElement(surfaceMember, _nsc('gml:Polygon'), attrib={_nsc('gml:id'): scene + '_7'})
        exterior = etree.SubElement(Polygon, _nsc('gml:exterior'))
        LinearRing = etree.SubElement(exterior, _nsc('gml:LinearRing'))
        posList = etree.SubElement(LinearRing, _nsc('gml:posList'), attrib={'uom': 'deg'})
        posList.text = meta['source'][uid]['geom_xml_envelop']
        
        centerOf = etree.SubElement(Footprint, _nsc('eop:centerOf'))
        Point = etree.SubElement(centerOf, _nsc('gml:Point'), attrib={_nsc('gml:id'): scene + '_8'})
        pos = etree.SubElement(Point, _nsc('gml:pos'), attrib={'uom': 'deg'})
        pos.text = meta['source'][uid]['geom_xml_center']
        
        ################################################################################################################
        result = etree.SubElement(root, _nsc('om:result'))
        EarthObservationResult = etree.SubElement(result, _nsc('eop:EarthObservationResult'),
                                                  attrib={_nsc('gml:id'): scene + '_9'})
        
        product = etree.SubElement(EarthObservationResult, _nsc('eop:product'))
        ProductInformation = etree.SubElement(product, _nsc('nrb:ProductInformation'))
        fileName = etree.SubElement(ProductInformation, _nsc('eop:fileName'))
        ServiceReference = etree.SubElement(fileName, _nsc('ows:ServiceReference'),
                                            attrib={_nsc('xlink:href'): scene})
        RequestMessage = etree.SubElement(ServiceReference, _nsc('ows:RequestMessage'))
        version = etree.SubElement(ProductInformation, _nsc('eop:size'))
        
        ################################################################################################################
        metaDataProperty = etree.SubElement(root, _nsc('eop:metaDataProperty'))
        EarthObservationMetaData = etree.SubElement(metaDataProperty, _nsc('nrb:EarthObservationMetaData'))
        
        identifier = etree.SubElement(EarthObservationMetaData, _nsc('eop:identifier'))
        identifier.text = scene
        doi = etree.SubElement(EarthObservationMetaData, _nsc('eop:doi'))
        doi.text = meta['source'][uid]['doi']
        status = etree.SubElement(EarthObservationMetaData, _nsc('eop:status'))
        status.text = meta['source'][uid]['status']
        acquisitionType = etree.SubElement(EarthObservationMetaData, _nsc('eop:acquisitionType'))
        acquisitionType.text = meta['source'][uid]['acquisitionType']
        productType = etree.SubElement(EarthObservationMetaData, _nsc('nrb:productType'),
                                       attrib={'codeSpace': 'urn:esa:eop:Sentinel1:class'})
        productType.text = meta['source'][uid]['productType']
        
        processing = etree.SubElement(EarthObservationMetaData, _nsc('nrb:processing'))
        ProcessingInformation = etree.SubElement(processing, _nsc('nrb:ProcessingInformation'))
        processingCenter = etree.SubElement(ProcessingInformation, _nsc('eop:processingCenter'),
                                            attrib={'codeSpace': 'urn:esa:eop:Sentinel1:facility'})
        processingCenter.text = meta['source'][uid]['processingCenter']
        processingDate = etree.SubElement(ProcessingInformation, _nsc('eop:processingDate'))
        processingDate.text = meta['source'][uid]['processingDate']
        processorName = etree.SubElement(ProcessingInformation, _nsc('eop:processorName'))
        processorName.text = meta['source'][uid]['processorName']
        processorVersion = etree.SubElement(ProcessingInformation, _nsc('eop:processorVersion'))
        processorVersion.text = meta['source'][uid]['processorVersion']
        processingLevel = etree.SubElement(ProcessingInformation, _nsc('eop:processingLevel'))
        processingLevel.text = meta['common']['processingLevel']
        processingMode = etree.SubElement(ProcessingInformation, _nsc('eop:processingMode'))
        processingMode.text = meta['source'][uid]['processingMode']
        orbitStateVector = etree.SubElement(ProcessingInformation, _nsc('nrb:orbitStateVector'),
                                            attrib={'access': meta['source'][uid]['orbitDataAccess']})
        orbitStateVector.text = meta['source'][uid]['orbitStateVector']
        for swath in meta['source'][uid]['swaths']:
            azimuthLookBandwidth = etree.SubElement(ProcessingInformation, _nsc('nrb:azimuthLookBandwidth'),
                                                    attrib={'uom': 'Hz', 'beam': swath})
            azimuthLookBandwidth.text = str(meta['source'][uid]['azimuthLookBandwidth'][swath])
        for swath in meta['source'][uid]['swaths']:
            rangeLookBandwidth = etree.SubElement(ProcessingInformation, _nsc('nrb:rangeLookBandwidth'),
                                                  attrib={'uom': 'Hz', 'beam': swath})
            rangeLookBandwidth.text = str(meta['source'][uid]['rangeLookBandwidth'][swath])
        
        lutApplied = etree.SubElement(ProcessingInformation, _nsc('nrb:lutApplied'))
        lutApplied.text = meta['source'][uid]['lutApplied']
        
        performance = etree.SubElement(EarthObservationMetaData, _nsc('nrb:performance'))
        PerformanceIndicators = etree.SubElement(performance, _nsc('nrb:PerformanceIndicators'))
        noiseEquivalentIntensity = etree.SubElement(PerformanceIndicators, _nsc('nrb:noiseEquivalentIntensity'),
                                                    attrib={'uom': 'dB',
                                                            'type': str(meta['source'][uid]['perfNoiseEquivalentIntensityType'])})
        for pol in meta['common']['polarisationChannels']:
            estimatesMin = etree.SubElement(noiseEquivalentIntensity, _nsc('nrb:estimates'),
                                            attrib={'pol': pol, 'type': 'minimum'})
            estimatesMin.text = str(meta['source'][uid]['perfEstimates'][pol]['minimum'])
            estimatesMax = etree.SubElement(noiseEquivalentIntensity, _nsc('nrb:estimates'),
                                            attrib={'pol': pol, 'type': 'maximum'})
            estimatesMax.text = str(meta['source'][uid]['perfEstimates'][pol]['maximum'])
            estimatesMean = etree.SubElement(noiseEquivalentIntensity, _nsc('nrb:estimates'),
                                             attrib={'pol': pol, 'type': 'mean'})
            estimatesMean.text = str(meta['source'][uid]['perfEstimates'][pol]['mean'])
        equivalentNumberOfLooks = etree.SubElement(PerformanceIndicators, _nsc('nrb:equivalentNumberOfLooks'))
        equivalentNumberOfLooks.text = str(meta['source'][uid]['perfEquivalentNumberOfLooks'])
        peakSideLobeRatio = etree.SubElement(PerformanceIndicators, _nsc('nrb:peakSideLobeRatio'),
                                             attrib={'uom': 'dB'})
        peakSideLobeRatio.text = str(meta['source'][uid]['perfPeakSideLobeRatio'])
        integratedSideLobeRatio = etree.SubElement(PerformanceIndicators, _nsc('nrb:integratedSideLobeRatio'),
                                                   attrib={'uom': 'dB'})
        integratedSideLobeRatio.text = str(meta['source'][uid]['perfIntegratedSideLobeRatio'])
        
        azimuthNumberOfLooks = etree.SubElement(EarthObservationMetaData, _nsc('nrb:azimuthNumberOfLooks'))
        azimuthNumberOfLooks.text = meta['source'][uid]['azimuthNumberOfLooks']
        rangeNumberOfLooks = etree.SubElement(EarthObservationMetaData, _nsc('nrb:rangeNumberOfLooks'))
        rangeNumberOfLooks.text = meta['source'][uid]['rangeNumberOfLooks']
        dataGeometry = etree.SubElement(EarthObservationMetaData, _nsc('nrb:dataGeometry'))
        dataGeometry.text = meta['source'][uid]['dataGeometry']
        for swath in meta['source'][uid]['swaths']:
            azimuthResolution = etree.SubElement(EarthObservationMetaData, _nsc('nrb:azimuthResolution'),
                                                 attrib={'uom': 'm', 'beam': swath})
            azimuthResolution.text = meta['source'][uid]['azimuthResolution'][swath]
        for swath in meta['source'][uid]['swaths']:
            rangeResolution = etree.SubElement(EarthObservationMetaData, _nsc('nrb:rangeResolution'),
                                               attrib={'uom': 'm', 'beam': swath})
            rangeResolution.text = meta['source'][uid]['rangeResolution'][swath]
        azimuthPixelSpacing = etree.SubElement(EarthObservationMetaData, _nsc('nrb:azimuthPixelSpacing'),
                                               attrib={'uom': 'm'})
        azimuthPixelSpacing.text = meta['source'][uid]['azimuthPixelSpacing']
        rangePixelSpacing = etree.SubElement(EarthObservationMetaData, _nsc('nrb:rangePixelSpacing'),
                                             attrib={'uom': 'm'})
        rangePixelSpacing.text = meta['source'][uid]['rangePixelSpacing']
        polCalMatrices = etree.SubElement(EarthObservationMetaData, _nsc('nrb:polCalMatrices'),
                                          attrib={_nsc('xlink:href'): str(meta['source'][uid]['polCalMatrices'])})
        meanFaradayRotationAngle = etree.SubElement(EarthObservationMetaData, _nsc('nrb:meanFaradayRotationAngle'),
                                                    attrib={'uom': 'deg'})
        meanFaradayRotationAngle.text = meta['source'][uid]['faradayMeanRotationAngle']
        referenceFaradayRotation = etree.SubElement(EarthObservationMetaData, _nsc('nrb:referenceFaradayRotation'),
                                                    attrib={_nsc('xlink:href'): str(meta['source'][uid]['faradayRotationReference'])})
        ionosphereIndicator = etree.SubElement(EarthObservationMetaData, _nsc('nrb:ionosphereIndicator'))
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
