import os
import re
from copy import deepcopy
from lxml import etree
from datetime import datetime, timezone
from spatialist import Raster
from spatialist.ancillary import finder
from statistics import mean
from S1_NRB.metadata.mapping import ASSET_MAP, NS_MAP
from S1_NRB.metadata.extract import get_header_size


def parse(meta, target, assets, exist_ok=False):
    """
    Wrapper for :func:`~S1_NRB.metadata.xml.source_xml` and :func:`~S1_NRB.metadata.xml.product_xml`.
    
    Parameters
    ----------
    meta: dict
        Metadata dictionary generated with :func:`~S1_NRB.metadata.extract.meta_dict`.
    target: str
        A path pointing to the root directory of a product scene.
    assets: list[str]
        List of paths to all GeoTIFF and VRT assets of the currently processed ARD product.
    exist_ok: bool
        Do not create files if they already exist?
    """
    key = "s1-{}".format(meta['prod']['productName-short'].lower())
    nsmap = deepcopy(NS_MAP)
    nsmap[key] = nsmap.pop('placeholder')
    src_url = nsmap[key].replace('spec', key.split('-')[1]).replace('role', 'source')
    prod_url = nsmap[key].replace('spec', key.split('-')[1]).replace('role', 'product')
    
    nsmap.update({key: src_url})
    source_xml(meta=meta, target=target, nsmap=nsmap, ard_ns=key, exist_ok=exist_ok)
    nsmap.update({key: prod_url})
    product_xml(meta=meta, target=target, assets=assets, nsmap=nsmap, ard_ns=key,
                exist_ok=exist_ok)


def source_xml(meta, target, nsmap, ard_ns, exist_ok=False):
    """
    Function to generate source-level metadata for an ARD product in `OGC 10-157r4` compliant XML format.
    
    Parameters
    ----------
    meta: dict
        Metadata dictionary generated with :func:`~S1_NRB.metadata.extract.meta_dict`
    target: str
        A path pointing to the root directory of a product scene.
    nsmap: dict
        Dictionary listing abbreviation (key) and URI (value) of all necessary XML namespaces.
    ard_ns: str
        Abbreviation of the ARD namespace. E.g., `s1-nrb` for the NRB ARD product.
    exist_ok: bool
        Do not create files if they already exist?
    """
    metadir = os.path.join(target, 'source')
    for uid in list(meta['source'].keys()):
        scene = os.path.basename(meta['source'][uid]['filename']).split('.')[0]
        outname = os.path.join(metadir, '{}.xml'.format(scene))
        if os.path.isfile(outname) and exist_ok:
            continue
        print(outname)
        timeStart = datetime.strftime(meta['source'][uid]['timeStart'], '%Y-%m-%dT%H:%M:%S.%f')
        timeStop = datetime.strftime(meta['source'][uid]['timeStop'], '%Y-%m-%dT%H:%M:%S.%f')
        
        root = etree.Element(_nsc('_:EarthObservation', nsmap, ard_ns=ard_ns), nsmap=nsmap,
                             attrib={_nsc('gml:id', nsmap): scene + '_1'})
        _om_time(root=root, nsmap=nsmap, scene_id=scene, time_start=timeStart, time_stop=timeStop)
        _om_procedure(root=root, nsmap=nsmap, ard_ns=ard_ns, scene_id=scene, meta=meta, uid=uid, prod=False)
        observedProperty = etree.SubElement(root, _nsc('om:observedProperty', nsmap),
                                            attrib={'nilReason': 'inapplicable'})
        _om_feature_of_interest(root=root, nsmap=nsmap, scene_id=scene,
                                extent=meta['source'][uid]['geom_xml_envelop'],
                                center=meta['source'][uid]['geom_xml_center'])
        
        ################################################################################################################
        result = etree.SubElement(root, _nsc('om:result', nsmap))
        earthObservationResult = etree.SubElement(result, _nsc('eop:EarthObservationResult', nsmap),
                                                  attrib={_nsc('gml:id', nsmap): scene + '_9'})
        product = etree.SubElement(earthObservationResult, _nsc('eop:product', nsmap))
        productInformation = etree.SubElement(product, _nsc('_:ProductInformation', nsmap, ard_ns=ard_ns))
        fileName = etree.SubElement(productInformation, _nsc('eop:fileName', nsmap))
        serviceReference = etree.SubElement(fileName, _nsc('ows:ServiceReference', nsmap),
                                            attrib={_nsc('xlink:href', nsmap): scene})
        requestMessage = etree.SubElement(serviceReference, _nsc('ows:RequestMessage', nsmap))
        
        org_src_files_dir = os.path.join(metadir, uid)
        if os.path.isdir(org_src_files_dir):
            org_src_files = finder(target=org_src_files_dir, matchlist=['*.safe', '*.xml'], foldermode=0)
            if len(org_src_files) > 0:
                for file in org_src_files:
                    href = './' + os.path.relpath(file, metadir).replace('\\', '/')
                    product = etree.SubElement(earthObservationResult, _nsc('eop:product', nsmap))
                    productInformation = etree.SubElement(product, _nsc('_:ProductInformation', nsmap, ard_ns=ard_ns))
                    fileName = etree.SubElement(productInformation, _nsc('eop:fileName', nsmap))
                    serviceReference = etree.SubElement(fileName, _nsc('ows:ServiceReference', nsmap),
                                                        attrib={_nsc('xlink:href', nsmap): href})
                    requestMessage = etree.SubElement(serviceReference, _nsc('ows:RequestMessage', nsmap))
                    dataFormat = etree.SubElement(productInformation, _nsc('_:dataFormat', nsmap, ard_ns=ard_ns))
                    dataFormat.text = 'XML'
        ################################################################################################################
        metaDataProperty = etree.SubElement(root, _nsc('eop:metaDataProperty', nsmap))
        earthObservationMetaData = etree.SubElement(metaDataProperty, _nsc('_:EarthObservationMetaData', nsmap,
                                                                           ard_ns=ard_ns))
        
        identifier = etree.SubElement(earthObservationMetaData, _nsc('eop:identifier', nsmap))
        identifier.text = scene
        doi = etree.SubElement(earthObservationMetaData, _nsc('eop:doi', nsmap))
        doi.text = meta['source'][uid]['doi']
        acquisitionType = etree.SubElement(earthObservationMetaData, _nsc('eop:acquisitionType', nsmap))
        acquisitionType.text = meta['source'][uid]['acquisitionType']
        status = etree.SubElement(earthObservationMetaData, _nsc('eop:status', nsmap))
        status.text = meta['source'][uid]['status']
        
        processing = etree.SubElement(earthObservationMetaData, _nsc('eop:processing', nsmap))
        processingInformation = etree.SubElement(processing, _nsc('_:ProcessingInformation', nsmap, ard_ns=ard_ns))
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
        processingLevel = etree.SubElement(processingInformation, _nsc('_:processingLevel', nsmap, ard_ns=ard_ns))
        processingLevel.text = meta['common']['processingLevel']
        orbitDataSource = etree.SubElement(processingInformation, _nsc('_:orbitDataSource', nsmap, ard_ns=ard_ns))
        orbitDataSource.text = meta['source'][uid]['orbitDataSource'].upper()
        orbitStateVector = etree.SubElement(processingInformation, _nsc('_:orbitStateVector', nsmap, ard_ns=ard_ns),
                                            attrib={'access': meta['source'][uid]['orbitDataAccess']})
        orbitStateVector.text = meta['source'][uid]['orbitStateVector']
        for swath in meta['source'][uid]['swaths']:
            azimuthLookBandwidth = etree.SubElement(processingInformation, _nsc('_:azimuthLookBandwidth', nsmap,
                                                                                ard_ns=ard_ns),
                                                    attrib={'uom': 'Hz', 'beam': swath})
            azimuthLookBandwidth.text = str(meta['source'][uid]['azimuthLookBandwidth'][swath])
        for swath in meta['source'][uid]['swaths']:  # removal will change order in output file
            rangeLookBandwidth = etree.SubElement(processingInformation, _nsc('_:rangeLookBandwidth', nsmap,
                                                                              ard_ns=ard_ns),
                                                  attrib={'uom': 'Hz', 'beam': swath})
            rangeLookBandwidth.text = str(meta['source'][uid]['rangeLookBandwidth'][swath])
        lutApplied = etree.SubElement(processingInformation, _nsc('_:lutApplied', nsmap, ard_ns=ard_ns))
        lutApplied.text = meta['source'][uid]['lutApplied']
        
        productType = etree.SubElement(earthObservationMetaData, _nsc('_:productType', nsmap, ard_ns=ard_ns),
                                       attrib={'codeSpace': 'urn:esa:eop:Sentinel1:class'})
        productType.text = meta['source'][uid]['productType']
        dataGeometry = etree.SubElement(earthObservationMetaData,
                                        _nsc('_:dataGeometry', nsmap, ard_ns=ard_ns))
        dataGeometry.text = meta['source'][uid]['dataGeometry']
        for swath in meta['source'][uid]['swaths']:
            azimuthNumberOfLooks = etree.SubElement(earthObservationMetaData,
                                                    _nsc('_:azimuthNumberOfLooks', nsmap, ard_ns=ard_ns),
                                                    attrib={'beam': swath})
            azimuthNumberOfLooks.text = str(meta['source'][uid]['azimuthNumberOfLooks'][swath])
        for swath in meta['source'][uid]['swaths']:  # removal will change order in output file
            rangeNumberOfLooks = etree.SubElement(earthObservationMetaData,
                                                  _nsc('_:rangeNumberOfLooks', nsmap, ard_ns=ard_ns),
                                                  attrib={'beam': swath})
            rangeNumberOfLooks.text = str(meta['source'][uid]['rangeNumberOfLooks'][swath])
        for swath in meta['source'][uid]['swaths']:
            azimuthResolution = etree.SubElement(earthObservationMetaData,
                                                 _nsc('_:azimuthResolution', nsmap, ard_ns=ard_ns),
                                                 attrib={'uom': 'm', 'beam': swath})
            azimuthResolution.text = str(meta['source'][uid]['azimuthResolution'][swath])
        for swath in meta['source'][uid]['swaths']:  # removal will change order in output file
            rangeResolution = etree.SubElement(earthObservationMetaData,
                                               _nsc('_:rangeResolution', nsmap, ard_ns=ard_ns),
                                               attrib={'uom': 'm', 'beam': swath})
            rangeResolution.text = str(meta['source'][uid]['rangeResolution'][swath])
        azimuthPixelSpacing = etree.SubElement(earthObservationMetaData, _nsc('_:azimuthPixelSpacing', nsmap,
                                                                              ard_ns=ard_ns),
                                               attrib={'uom': 'm'})
        azimuthPixelSpacing.text = str(mean(meta['source'][uid]['azimuthPixelSpacing'].values()))
        rangePixelSpacing = etree.SubElement(earthObservationMetaData, _nsc('_:rangePixelSpacing', nsmap,
                                                                            ard_ns=ard_ns),
                                             attrib={'uom': 'm'})
        rangePixelSpacing.text = str(mean(meta['source'][uid]['rangePixelSpacing'].values()))
        
        performance = etree.SubElement(earthObservationMetaData, _nsc('_:performance', nsmap, ard_ns=ard_ns))
        performanceIndicators = etree.SubElement(performance, _nsc('_:PerformanceIndicators', nsmap, ard_ns=ard_ns))
        noiseEquivalentIntensityType = etree.SubElement(performanceIndicators,
                                                        _nsc('_:noiseEquivalentIntensityType', nsmap, ard_ns=ard_ns),
                                                        attrib={'uom': 'dB'})
        noiseEquivalentIntensityType.text = str(meta['source'][uid]['perfNoiseEquivalentIntensityType'])
        for pol in meta['common']['polarisationChannels']:
            estimatesMin = etree.SubElement(performanceIndicators, _nsc('_:estimates', nsmap, ard_ns=ard_ns),
                                            attrib={'pol': pol, 'type': 'minimum'})
            estimatesMin.text = str(meta['source'][uid]['perfEstimates'][pol]['minimum'])
            estimatesMax = etree.SubElement(performanceIndicators, _nsc('_:estimates', nsmap, ard_ns=ard_ns),
                                            attrib={'pol': pol, 'type': 'maximum'})
            estimatesMax.text = str(meta['source'][uid]['perfEstimates'][pol]['maximum'])
            estimatesMean = etree.SubElement(performanceIndicators, _nsc('_:estimates', nsmap, ard_ns=ard_ns),
                                             attrib={'pol': pol, 'type': 'mean'})
            estimatesMean.text = str(meta['source'][uid]['perfEstimates'][pol]['mean'])
        equivalentNumberOfLooks = etree.SubElement(performanceIndicators, _nsc('_:equivalentNumberOfLooks', nsmap,
                                                                               ard_ns=ard_ns))
        equivalentNumberOfLooks.text = str(meta['source'][uid]['perfEquivalentNumberOfLooks'])
        peakSideLobeRatio = etree.SubElement(performanceIndicators, _nsc('_:peakSideLobeRatio', nsmap, ard_ns=ard_ns),
                                             attrib={'uom': 'dB'})
        peakSideLobeRatio.text = str(meta['source'][uid]['perfPeakSideLobeRatio'])
        integratedSideLobeRatio = etree.SubElement(performanceIndicators, _nsc('_:integratedSideLobeRatio', nsmap,
                                                                               ard_ns=ard_ns),
                                                   attrib={'uom': 'dB'})
        integratedSideLobeRatio.text = str(meta['source'][uid]['perfIntegratedSideLobeRatio'])
        
        polCalMatrices = etree.SubElement(earthObservationMetaData, _nsc('_:polCalMatrices', nsmap, ard_ns=ard_ns),
                                          attrib={
                                              _nsc('xlink:href', nsmap): str(meta['source'][uid]['polCalMatrices'])})
        meanFaradayRotationAngle = etree.SubElement(earthObservationMetaData,
                                                    _nsc('_:meanFaradayRotationAngle', nsmap, ard_ns=ard_ns),
                                                    attrib={'uom': 'deg'})
        meanFaradayRotationAngle.text = meta['source'][uid]['faradayMeanRotationAngle']
        faraday_ref = str(meta['source'][uid]['faradayRotationReference'])
        referenceFaradayRotation = etree.SubElement(earthObservationMetaData,
                                                    _nsc('_:referenceFaradayRotation', nsmap, ard_ns=ard_ns),
                                                    attrib={_nsc('xlink:href', nsmap): faraday_ref})
        ionosphereIndicator = etree.SubElement(earthObservationMetaData, _nsc('_:ionosphereIndicator', nsmap,
                                                                              ard_ns=ard_ns))
        ionosphereIndicator.text = meta['source'][uid]['ionosphereIndicator']
        
        ################################################################################################################
        etree.indent(root)
        tree = etree.ElementTree(root)
        tree.write(outname, pretty_print=True, xml_declaration=True, encoding='utf-8')


def product_xml(meta, target, assets, nsmap, ard_ns, exist_ok=False):
    """
    Function to generate product-level metadata for an ARD product in `OGC 10-157r4` compliant XML format.
    
    Parameters
    ----------
    meta: dict
        Metadata dictionary generated with :func:`~S1_NRB.metadata.extract.meta_dict`
    target: str
        A path pointing to the root directory of a product scene.
    assets: list[str]
        List of paths to all GeoTIFF and VRT assets of the currently processed ARD product.
    nsmap: dict
        Dictionary listing abbreviation (key) and URI (value) of all necessary XML namespaces.
    ard_ns: str
        Abbreviation of the ARD namespace. E.g., `s1-nrb` for the NRB ARD product.
    exist_ok: bool
        Do not create files if they already exist?
    """
    scene_id = os.path.basename(target)
    outname = os.path.join(target, '{}.xml'.format(scene_id))
    if os.path.isfile(outname) and exist_ok:
        return
    print(outname)
    timeCreated = datetime.strftime(meta['prod']['timeCreated'], '%Y-%m-%dT%H:%M:%S.%f')
    timeStart = datetime.strftime(meta['prod']['timeStart'], '%Y-%m-%dT%H:%M:%S.%f')
    timeStop = datetime.strftime(meta['prod']['timeStop'], '%Y-%m-%dT%H:%M:%S.%f')
    
    root = etree.Element(_nsc('_:EarthObservation', nsmap, ard_ns=ard_ns), nsmap=nsmap,
                         attrib={_nsc('gml:id', nsmap): scene_id + '_1'})
    _om_time(root=root, nsmap=nsmap, scene_id=scene_id, time_start=timeStart, time_stop=timeStop)
    _om_procedure(root=root, nsmap=nsmap, ard_ns=ard_ns, scene_id=scene_id, meta=meta, prod=True)
    observedProperty = etree.SubElement(root, _nsc('om:observedProperty', nsmap),
                                        attrib={'nilReason': 'inapplicable'})
    _om_feature_of_interest(root=root, nsmap=nsmap, scene_id=scene_id,
                            extent=meta['prod']['geom_xml_envelope'],
                            center=meta['prod']['geom_xml_center'])
    
    ####################################################################################################################
    result = etree.SubElement(root, _nsc('om:result', nsmap))
    earthObservationResult = etree.SubElement(result, _nsc('eop:EarthObservationResult', nsmap),
                                              attrib={_nsc('gml:id', nsmap): scene_id + '_9'})
    product = etree.SubElement(earthObservationResult, _nsc('eop:product', nsmap))
    productInformation = etree.SubElement(product, _nsc('_:ProductInformation', nsmap, ard_ns=ard_ns))
    fileName = etree.SubElement(productInformation, _nsc('eop:fileName', nsmap))
    serviceReference = etree.SubElement(fileName, _nsc('ows:ServiceReference', nsmap),
                                        attrib={_nsc('xlink:href', nsmap): scene_id})
    requestMessage = etree.SubElement(serviceReference, _nsc('ows:RequestMessage', nsmap))
    
    for asset in assets:
        relpath = './' + os.path.relpath(asset, target).replace('\\', '/')
        
        no_data = None
        header_size = None
        data_format = 'VRT'
        byte_order = None
        data_type = None
        z_error = None
        if asset.endswith('.tif'):
            with Raster(asset) as ras:
                no_data = str(ras.nodata)
            header_size = str(get_header_size(tif=asset))
            data_format = 'COG'
            byte_order = 'little-endian'
            data_type = 'FLOAT 32'
            prefix = '[0-9a-z]{5}-'
            match = re.search(prefix + f"({'|'.join(meta['prod']['compression_zerrors'].keys())})",
                              os.path.basename(asset))
            if match is not None:
                k = match.group()
                k = k.removeprefix(re.search(prefix, k).group())
                z_error = str(meta['prod']['compression_zerrors'][k])
        
        product = etree.SubElement(earthObservationResult, _nsc('eop:product', nsmap))
        productInformation = etree.SubElement(product, _nsc('_:ProductInformation', nsmap, ard_ns=ard_ns))
        fileName = etree.SubElement(productInformation, _nsc('eop:fileName', nsmap))
        serviceReference = etree.SubElement(fileName, _nsc('ows:ServiceReference', nsmap),
                                            attrib={_nsc('xlink:href', nsmap): relpath})
        requestMessage = etree.SubElement(serviceReference, _nsc('ows:RequestMessage', nsmap))
        
        size = etree.SubElement(productInformation, _nsc('eop:size', nsmap), attrib={'uom': 'bytes'})
        size.text = str(os.path.getsize(asset))
        if header_size is not None:
            headerSize = etree.SubElement(productInformation, _nsc('_:headerSize', nsmap, ard_ns=ard_ns),
                                          attrib={'uom': 'bytes'})
            headerSize.text = header_size
        if byte_order is not None:
            byteOrder = etree.SubElement(productInformation, _nsc('_:byteOrder', nsmap, ard_ns=ard_ns))
            byteOrder.text = byte_order
        dataFormat = etree.SubElement(productInformation, _nsc('_:dataFormat', nsmap, ard_ns=ard_ns))
        dataFormat.text = data_format
        if data_type is not None:
            dataType = etree.SubElement(productInformation, _nsc('_:dataType', nsmap, ard_ns=ard_ns))
            dataType.text = data_type.split()[0]
            bitsPerSample = etree.SubElement(productInformation, _nsc('_:bitsPerSample', nsmap, ard_ns=ard_ns))
            bitsPerSample.text = data_type.split()[1]
        if no_data is not None:
            noDataVal = etree.SubElement(productInformation, _nsc('_:noDataValue', nsmap, ard_ns=ard_ns))
            noDataVal.text = no_data
        if z_error is not None:
            compressionType = etree.SubElement(productInformation, _nsc('_:compressionType', nsmap, ard_ns=ard_ns))
            compressionType.text = meta['prod']['compression_type']
            compressionzError = etree.SubElement(productInformation, _nsc('_:compressionZError', nsmap, ard_ns=ard_ns))
            compressionzError.text = z_error
        
        if 'annotation' in asset:
            key = re.search('-[a-z]{2}(?:-[a-z]{2}|).tif', asset).group()
            np_pat = '-np-[vh]{2}.tif'
            if re.search(np_pat, key) is not None:
                key = np_pat
            
            sampleType = etree.SubElement(productInformation, _nsc('_:sampleType', nsmap, ard_ns=ard_ns),
                                          attrib={'uom': 'unitless' if ASSET_MAP[key]['unit'] is None else
                                          ASSET_MAP[key]['unit']})
            sampleType.text = ASSET_MAP[key]['type']
            
            if key in ['-dm.tif', '-id.tif']:
                dataType.text = 'UINT'
                bitsPerSample.text = '8'
                
                if key == '-dm.tif':
                    with Raster(asset) as dm_ras:
                        band_descr = [dm_ras.raster.GetRasterBand(band).GetDescription() for band in
                                      range(1, dm_ras.bands + 1)]
                    samples = [x for x in band_descr if x in ASSET_MAP[key]['allowed']]
                    for i, sample in enumerate(samples):
                        bitValue = etree.SubElement(productInformation, _nsc('_:bitValue', nsmap, ard_ns=ard_ns),
                                                    attrib={'band': str(i + 1),
                                                            'name': sample})
                        bitValue.text = '1'
                else:  # key == '-id.tif'
                    src_list = list(meta['source'].keys())
                    src_target = [os.path.basename(meta['source'][src]['filename']).replace('.SAFE',
                                                                                            '').replace('.zip', '')
                                  for src in src_list]
                    for i, s in enumerate(src_target):
                        bitValue = etree.SubElement(productInformation, _nsc('_:bitValue', nsmap, ard_ns=ard_ns),
                                                    attrib={'band': '1', 'name': s})
                        bitValue.text = str(i + 1)
            
            if key == '-ei.tif':
                ellipsoidalHeight = etree.SubElement(productInformation, _nsc('_:ellipsoidalHeight', nsmap,
                                                                              ard_ns=ard_ns),
                                                     attrib={'uom': 'm'})
                ellipsoidalHeight.text = meta['prod']['ellipsoidalHeight']
        
        if 'measurement' in asset and not asset.endswith('.vrt'):
            creationTime = etree.SubElement(productInformation, _nsc('_:creationTime', nsmap, ard_ns=ard_ns))
            creationTime.text = datetime.fromtimestamp(os.path.getctime(asset), tz=timezone.utc).isoformat()
            polarization = etree.SubElement(productInformation, _nsc('_:polarization', nsmap, ard_ns=ard_ns))
            polarization.text = re.search('-[vh]{2}', relpath).group().removeprefix('-').upper()
            numBorderPixels = etree.SubElement(productInformation, _nsc('_:numBorderPixels', nsmap, ard_ns=ard_ns))
            numBorderPixels.text = str(meta['prod']['numBorderPixels'])
    
    ####################################################################################################################
    metaDataProperty = etree.SubElement(root, _nsc('eop:metaDataProperty', nsmap))
    earthObservationMetaData = etree.SubElement(metaDataProperty, _nsc('_:EarthObservationMetaData', nsmap,
                                                                       ard_ns=ard_ns))
    
    identifier = etree.SubElement(earthObservationMetaData, _nsc('eop:identifier', nsmap))
    identifier.text = scene_id
    if meta['prod']['doi'] is not None:
        doi = etree.SubElement(earthObservationMetaData, _nsc('eop:doi', nsmap))
        doi.text = meta['prod']['doi']
    acquisitionType = etree.SubElement(earthObservationMetaData, _nsc('eop:acquisitionType', nsmap))
    acquisitionType.text = meta['prod']['acquisitionType']
    status = etree.SubElement(earthObservationMetaData, _nsc('eop:status', nsmap))
    status.text = meta['prod']['status']
    
    processing = etree.SubElement(earthObservationMetaData, _nsc('eop:processing', nsmap))
    processingInformation = etree.SubElement(processing, _nsc('_:ProcessingInformation', nsmap, ard_ns=ard_ns))
    if meta['prod']['processingCenter'] is not None:
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
    processingLevel = etree.SubElement(processingInformation, _nsc('_:processingLevel', nsmap, ard_ns=ard_ns))
    processingLevel.text = meta['common']['processingLevel']
    for src in list(meta['source'].keys()):
        src_path = '{}.xml'.format(os.path.basename(meta['source'][src]['filename']).split('.')[0])
        src_target = os.path.join('./source', src_path).replace('\\', '/')
        sourceProduct = etree.SubElement(processingInformation, _nsc('_:sourceProduct', nsmap, ard_ns=ard_ns),
                                         attrib={_nsc('xlink:href', nsmap): src_target})
    auxData1 = etree.SubElement(processingInformation, _nsc('_:auxiliaryDataSetFileName', nsmap, ard_ns=ard_ns),
                                attrib={_nsc('xlink:href', nsmap): meta['prod']['ancillaryData_KML']})
    speckleFilterApplied = etree.SubElement(processingInformation, _nsc('_:speckleFilterApplied', nsmap,
                                                                        ard_ns=ard_ns))
    speckleFilterApplied.text = str(meta['prod']['speckleFilterApplied']).lower()
    noiseRemovalApplied = etree.SubElement(processingInformation, _nsc('_:noiseRemovalApplied', nsmap, ard_ns=ard_ns))
    noiseRemovalApplied.text = str(meta['prod']['noiseRemovalApplied']).lower()
    if meta['prod']['noiseRemovalApplied']:
        noiseRemovalAlgorithm = etree.SubElement(processingInformation,
                                                 _nsc('_:noiseRemovalAlgorithm', nsmap, ard_ns=ard_ns),
                                                 attrib={_nsc('xlink:href', nsmap):
                                                             meta['prod']['noiseRemovalAlgorithm']})
    if meta['prod']['RTCAlgorithm'] is not None:
        rtcAlgorithm = etree.SubElement(processingInformation, _nsc('_:RTCAlgorithm', nsmap, ard_ns=ard_ns),
                                        attrib={_nsc('xlink:href', nsmap): meta['prod']['RTCAlgorithm']})
    if meta['prod']['windNormBackscatterMeasurement'] is not None:
        windNormBackscatterMeasurement = etree.SubElement(processingInformation,
                                                          _nsc('_:windNormBackscatterMeasurement',
                                                               nsmap, ard_ns=ard_ns))
        windNormBackscatterMeasurement.text = meta['prod']['windNormBackscatterMeasurement']
        windNormBackscatterConvention = etree.SubElement(processingInformation,
                                                         _nsc('_:windNormBackscatterConvention', nsmap, ard_ns=ard_ns))
        windNormBackscatterConvention.text = meta['prod']['windNormBackscatterConvention']
        windNormReferenceDirection = etree.SubElement(processingInformation,
                                                      _nsc('_:windNormReferenceDirection', nsmap, ard_ns=ard_ns),
                                                      attrib={'uom': 'deg'})
        windNormReferenceDirection.text = str(meta['prod']['windNormReferenceDirection'])
        
        windNormReferenceModel = etree.SubElement(processingInformation, _nsc('_:windNormReferenceModel', nsmap,
                                                                              ard_ns=ard_ns),
                                                  attrib={_nsc('xlink:href', nsmap):
                                                              meta['prod']['windNormReferenceModel']})
        windNormReferenceSpeed = etree.SubElement(processingInformation,
                                                  _nsc('_:windNormReferenceSpeed', nsmap, ard_ns=ard_ns),
                                                  attrib={'uom': 'm_s'})
        windNormReferenceSpeed.text = str(meta['prod']['windNormReferenceSpeed'])
        windNormReferenceType = etree.SubElement(processingInformation,
                                                 _nsc('_:windNormReferenceType', nsmap, ard_ns=ard_ns))
        windNormReferenceType.text = meta['prod']['windNormReferenceType']
    geoCorrAlgorithm = etree.SubElement(processingInformation, _nsc('_:geoCorrAlgorithm', nsmap, ard_ns=ard_ns),
                                        attrib={_nsc('xlink:href', nsmap): meta['prod']['geoCorrAlgorithm']})
    geoCorrResamplingMethod = etree.SubElement(processingInformation, _nsc('_:geoCorrResamplingAlgorithm', nsmap,
                                                                           ard_ns=ard_ns))
    geoCorrResamplingMethod.text = meta['prod']['geoCorrResamplingMethod'].upper()
    demReference = etree.SubElement(processingInformation, _nsc('_:DEMReference', nsmap, ard_ns=ard_ns),
                                    attrib={'name': meta['prod']['demName'],
                                            'dem': meta['prod']['demType'],
                                            _nsc('xlink:href', nsmap): meta['prod']['demReference']})
    demResamplingMethod = etree.SubElement(processingInformation, _nsc('_:DEMResamplingMethod', nsmap, ard_ns=ard_ns))
    demResamplingMethod.text = meta['prod']['demResamplingMethod'].upper()
    demAccess = etree.SubElement(processingInformation, _nsc('_:DEMAccess', nsmap, ard_ns=ard_ns),
                                 attrib={_nsc('xlink:href', nsmap): meta['prod']['demAccess']})
    demGSD = etree.SubElement(processingInformation, _nsc('_:DEMGroundSamplingDistance', nsmap, ard_ns=ard_ns),
                              attrib={'uom': meta['prod']['demGSD'].split()[1]})
    demGSD.text = meta['prod']['demGSD'].split()[0]
    egmReference = etree.SubElement(processingInformation, _nsc('_:EGMReference', nsmap, ard_ns=ard_ns),
                                    attrib={_nsc('xlink:href', nsmap): meta['prod']['demEGMReference']})
    egmResamplingMethod = etree.SubElement(processingInformation, _nsc('_:EGMResamplingMethod', nsmap,
                                                                       ard_ns=ard_ns))
    egmResamplingMethod.text = meta['prod']['demEGMResamplingMethod'].upper()
    
    productType = etree.SubElement(earthObservationMetaData, _nsc('_:productType', nsmap, ard_ns=ard_ns),
                                   attrib={'codeSpace': 'urn:esa:eop:Sentinel1:class'})
    productType.text = meta['prod']['productName-short']
    refDoc = etree.SubElement(earthObservationMetaData, _nsc('_:refDoc', nsmap, ard_ns=ard_ns),
                              attrib={'name': meta['prod']['productName'],
                                      'version': meta['prod']['card4l-version'],
                                      _nsc('xlink:href', nsmap): meta['prod']['card4l-link']})
    azimuthNumberOfLooks = etree.SubElement(earthObservationMetaData, _nsc('_:azimuthNumberOfLooks', nsmap,
                                                                           ard_ns=ard_ns))
    azimuthNumberOfLooks.text = str(meta['prod']['azimuthNumberOfLooks'])
    rangeNumberOfLooks = etree.SubElement(earthObservationMetaData, _nsc('_:rangeNumberOfLooks', nsmap,
                                                                         ard_ns=ard_ns))
    rangeNumberOfLooks.text = str(meta['prod']['rangeNumberOfLooks'])
    equivalentNumberLooks = etree.SubElement(earthObservationMetaData, _nsc('_:equivalentNumberOfLooks', nsmap,
                                                                            ard_ns=ard_ns))
    equivalentNumberLooks.text = str(meta['prod']['equivalentNumberLooks'])
    radiometricAccuracyRelative = etree.SubElement(earthObservationMetaData,
                                                   _nsc('_:radiometricAccuracyRelative', nsmap, ard_ns=ard_ns),
                                                   attrib={'uom': 'dB'})
    radiometricAccuracyRelative.text = meta['prod']['radiometricAccuracyRelative']
    radiometricAccuracyAbsolute = etree.SubElement(earthObservationMetaData,
                                                   _nsc('_:radiometricAccuracyAbsolute', nsmap, ard_ns=ard_ns),
                                                   attrib={'uom': 'dB'})
    radiometricAccuracyAbsolute.text = meta['prod']['radiometricAccuracyAbsolute']
    radacc_ref = str(meta['prod']['radiometricAccuracyReference'])
    radiometricAccuracyReference = etree.SubElement(earthObservationMetaData,
                                                    _nsc('_:radiometricAccuracyReference', nsmap, ard_ns=ard_ns),
                                                    attrib={_nsc('xlink:href', nsmap): radacc_ref})
    geoCorrAccuracyType = etree.SubElement(earthObservationMetaData, _nsc('_:geoCorrAccuracyType', nsmap,
                                                                          ard_ns=ard_ns))
    geoCorrAccuracyType.text = meta['prod']['geoCorrAccuracyType']
    geoCorrAccuracyNorthernSTDev = etree.SubElement(earthObservationMetaData,
                                                    _nsc('_:geoCorrAccuracyNorthernSTDev', nsmap, ard_ns=ard_ns),
                                                    attrib={'uom': 'm'})
    geoCorrAccuracyNorthernSTDev.text = meta['prod']['geoCorrAccuracyNorthernSTDev']
    geoCorrAccuracyEasternSTDev = etree.SubElement(earthObservationMetaData,
                                                   _nsc('_:geoCorrAccuracyEasternSTDev', nsmap, ard_ns=ard_ns),
                                                   attrib={'uom': 'm'})
    geoCorrAccuracyEasternSTDev.text = meta['prod']['geoCorrAccuracyEasternSTDev']
    geoCorrAccuracyNorthernBias = etree.SubElement(earthObservationMetaData,
                                                   _nsc('_:geoCorrAccuracyNorthernBias', nsmap, ard_ns=ard_ns),
                                                   attrib={'uom': 'm'})
    geoCorrAccuracyNorthernBias.text = meta['prod']['geoCorrAccuracyNorthernBias']
    geoCorrAccuracyEasternBias = etree.SubElement(earthObservationMetaData,
                                                  _nsc('_:geoCorrAccuracyEasternBias', nsmap, ard_ns=ard_ns),
                                                  attrib={'uom': 'm'})
    geoCorrAccuracyEasternBias.text = meta['prod']['geoCorrAccuracyEasternBias']
    geoCorrAccuracy_rRMSE = etree.SubElement(earthObservationMetaData,
                                             _nsc('_:geoCorrAccuracy_rRMSE', nsmap, ard_ns=ard_ns),
                                             attrib={'uom': 'm'})
    geoCorrAccuracy_rRMSE.text = str(meta['prod']['geoCorrAccuracy_rRMSE'])
    geoacc_ref = meta['prod']['geoCorrAccuracyReference']
    if geoacc_ref is not None:
        geoCorrAccuracyReference = etree.SubElement(earthObservationMetaData,
                                                    _nsc('_:geoCorrAccuracyReference', nsmap, ard_ns=ard_ns),
                                                    attrib={_nsc('xlink:href', nsmap): geoacc_ref})
    numLines = etree.SubElement(earthObservationMetaData, _nsc('_:numLines', nsmap, ard_ns=ard_ns))
    numLines.text = meta['prod']['numLines']
    numPixelsPerLine = etree.SubElement(earthObservationMetaData, _nsc('_:numPixelsPerLine', nsmap, ard_ns=ard_ns))
    numPixelsPerLine.text = meta['prod']['numPixelsPerLine']
    columnSpacing = etree.SubElement(earthObservationMetaData, _nsc('_:columnSpacing', nsmap, ard_ns=ard_ns),
                                     attrib={'uom': 'm'})
    columnSpacing.text = meta['prod']['pxSpacingColumn']
    rowSpacing = etree.SubElement(earthObservationMetaData, _nsc('_:rowSpacing', nsmap, ard_ns=ard_ns),
                                  attrib={'uom': 'm'})
    rowSpacing.text = meta['prod']['pxSpacingRow']
    pixelCoordinateConvention = etree.SubElement(earthObservationMetaData,
                                                 _nsc('_:pixelCoordinateConvention', nsmap, ard_ns=ard_ns))
    pixelCoordinateConvention.text = meta['prod']['pixelCoordinateConvention']
    backscatterMeasurement = etree.SubElement(earthObservationMetaData, _nsc('_:backscatterMeasurement', nsmap,
                                                                             ard_ns=ard_ns))
    backscatterMeasurement.text = meta['prod']['backscatterMeasurement']
    backscatterConvention = etree.SubElement(earthObservationMetaData, _nsc('_:backscatterConvention', nsmap,
                                                                            ard_ns=ard_ns))
    backscatterConvention.text = meta['prod']['backscatterConvention']
    backscatterConversionEq = etree.SubElement(earthObservationMetaData, _nsc('_:backscatterConversionEq', nsmap,
                                                                              ard_ns=ard_ns),
                                               attrib={'uom': 'dB'})
    backscatterConversionEq.text = meta['prod']['backscatterConversionEq']
    griddingConvention = etree.SubElement(earthObservationMetaData, _nsc('_:griddingConvention', nsmap, ard_ns=ard_ns),
                                          attrib={_nsc('xlink:href', nsmap): meta['prod']['griddingConventionURL']})
    mgrsID = etree.SubElement(earthObservationMetaData, _nsc('_:mgrsID', nsmap, ard_ns=ard_ns))
    mgrsID.text = meta['prod']['mgrsID']
    crsEPSG = etree.SubElement(earthObservationMetaData, _nsc('_:crsEPSG', nsmap, ard_ns=ard_ns),
                               attrib={'codeSpace': 'urn:esa:eop:crs'})
    crsEPSG.text = meta['prod']['crsEPSG']
    crsWKT = etree.SubElement(earthObservationMetaData, _nsc('_:crsWKT', nsmap, ard_ns=ard_ns))
    crsWKT.text = meta['prod']['crsWKT']
    
    ####################################################################################################################
    etree.indent(root)
    tree = etree.ElementTree(root)
    tree.write(outname, pretty_print=True, xml_declaration=True, encoding='utf-8')


def _nsc(text, nsmap, ard_ns=None):
    ns, key = text.split(':')
    if ard_ns is not None and ns == '_':
        ns = ard_ns
    return '{{{0}}}{1}'.format(nsmap[ns], key)


def _om_time(root, nsmap, scene_id, time_start, time_stop):
    """
    Creates the `om:phenomenonTime` and `om:resultTime` XML elements.
    
    Parameters
    ----------
    root: lxml.etree.Element
        Root XML element.
    nsmap: dict
        Dictionary listing abbreviation (key) and URI (value) of all necessary XML namespaces.
    scene_id: str
        Scene basename.
    time_start: str
        Start time of the scene acquisition.
    time_stop: str
        Stop time of the acquisition.
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


def _om_procedure(root, nsmap, ard_ns, scene_id, meta, uid=None, prod=True):
    """
    Creates the `om:procedure/eop:EarthObservationEquipment` XML elements and all relevant subelements for source and
    product metadata. Differences between source and product are controlled using the `prod=[True|False]` switch.
    
    Parameters
    ----------
    root: lxml.etree.Element
        Root XML element.
    nsmap: dict
        Dictionary listing abbreviation (key) and URI (value) of all necessary XML namespaces.
    ard_ns: str
        Abbreviation of the ARD namespace. E.g., `s1-nrb` for the NRB ARD product.
    scene_id: str
        Scene basename.
    meta: dict
        Metadata dictionary generated with :func:`~S1_NRB.metadata.extract.meta_dict`
    uid: str or None
        Unique identifier of a source SLC scene.
    prod: bool
        Return XML subelements for further usage in :func:`~S1_NRB.metadata.xml.product_xml` parsing function?
        Default is True. If False, the XML subelements for further usage in the :func:`~S1_NRB.metadata.xml.source_xml`
        parsing function will be returned.
    """
    procedure = etree.SubElement(root, _nsc('om:procedure', nsmap))
    earthObservationEquipment = etree.SubElement(procedure, _nsc('eop:EarthObservationEquipment', nsmap),
                                                 attrib={_nsc('gml:id', nsmap): scene_id + '_4'})
    
    # eop:platform
    platform0 = etree.SubElement(earthObservationEquipment, _nsc('eop:platform', nsmap))
    if prod:
        platform1 = etree.SubElement(platform0, _nsc('eop:Platform', nsmap))
    else:
        platform1 = etree.SubElement(platform0, _nsc('_:Platform', nsmap, ard_ns=ard_ns))
    shortName = etree.SubElement(platform1, _nsc('eop:shortName', nsmap))
    shortName.text = meta['common']['platformShortName'].upper()
    serialIdentifier = etree.SubElement(platform1, _nsc('eop:serialIdentifier', nsmap))
    serialIdentifier.text = meta['common']['platformIdentifier']
    if not prod:
        satReference = etree.SubElement(platform1, _nsc('_:satelliteReference', nsmap, ard_ns=ard_ns),
                                        attrib={_nsc('xlink:href', nsmap): meta['common']['platformReference']})
    
    # eop:instrument
    instrument0 = etree.SubElement(earthObservationEquipment, _nsc('eop:instrument', nsmap))
    instrument1 = etree.SubElement(instrument0, _nsc('eop:Instrument', nsmap))
    shortName = etree.SubElement(instrument1, _nsc('eop:shortName', nsmap))
    shortName.text = meta['common']['instrumentShortName']
    
    # eop:sensor
    sensor0 = etree.SubElement(earthObservationEquipment, _nsc('eop:sensor', nsmap))
    sensor1 = etree.SubElement(sensor0, _nsc('_:Sensor', nsmap, ard_ns=ard_ns))
    sensorType = etree.SubElement(sensor1, _nsc('eop:sensorType', nsmap))
    sensorType.text = meta['common']['sensorType']
    operationalMode = etree.SubElement(sensor1, _nsc('eop:operationalMode', nsmap),
                                       attrib={'codeSpace': 'urn:esa:eop:C-SAR:operationalMode'})
    operationalMode.text = meta['common']['operationalMode']
    swathIdentifier = etree.SubElement(sensor1, _nsc('eop:swathIdentifier', nsmap),
                                       attrib={'codeSpace': 'urn:esa:eop:C-SAR:swathIdentifier'})
    swathIdentifier.text = meta['common']['swathIdentifier']
    radarBand = etree.SubElement(sensor1, _nsc('_:radarBand', nsmap, ard_ns=ard_ns))
    radarBand.text = meta['common']['radarBand']
    if not prod:
        radarCenterFreq = etree.SubElement(sensor1, _nsc('_:radarCenterFrequency', nsmap, ard_ns=ard_ns),
                                           attrib={'uom': 'Hz'})
        radarCenterFreq.text = '{:.3e}'.format(meta['common']['radarCenterFreq'])
        sensorCalibration = etree.SubElement(sensor1, _nsc('_:sensorCalibration', nsmap, ard_ns=ard_ns),
                                             attrib={
                                                 _nsc('xlink:href', nsmap): meta['source'][uid]['sensorCalibration']})
    
    # eop:acquisitionParameters
    acquisitionParameters = etree.SubElement(earthObservationEquipment, _nsc('eop:acquisitionParameters', nsmap))
    acquisition = etree.SubElement(acquisitionParameters, _nsc('_:Acquisition', nsmap, ard_ns=ard_ns))
    orbitNumber = etree.SubElement(acquisition, _nsc('eop:orbitNumber', nsmap))
    orbitNumber.text = str(meta['common']['orbitNumber_abs'])
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
        instrumentAzimuthAngle.text = str(meta['source'][uid]['instrumentAzimuthAngle'])
    polarisationMode = etree.SubElement(acquisition, _nsc('sar:polarisationMode', nsmap))
    polarisationMode.text = meta['common']['polarisationMode']
    polarisationChannels = etree.SubElement(acquisition, _nsc('sar:polarisationChannels', nsmap))
    polarisationChannels.text = ', '.join(meta['common']['polarisationChannels'])
    if prod:
        numberOfAcquisitions = etree.SubElement(acquisition, _nsc('_:numberOfAcquisitions', nsmap, ard_ns=ard_ns))
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
        orbitMeanAltitude = etree.SubElement(acquisition, _nsc('_:orbitMeanAltitude', nsmap, ard_ns=ard_ns),
                                             attrib={'uom': 'm'})
        orbitMeanAltitude.text = meta['common']['orbitMeanAltitude']
        dataTakeID = etree.SubElement(acquisition, _nsc('_:dataTakeID', nsmap, ard_ns=ard_ns))
        dataTakeID.text = meta['source'][uid]['datatakeID']
        majorCycleID = etree.SubElement(acquisition, _nsc('_:majorCycleID', nsmap, ard_ns=ard_ns))
        majorCycleID.text = meta['source'][uid]['majorCycleID']


def _om_feature_of_interest(root, nsmap, scene_id, extent, center):
    """
    Creates the `om:featureOfInterest` XML elements.
    
    Parameters
    ----------
    root: lxml.etree.Element
        Root XML element.
    nsmap: dict
        Dictionary listing abbreviation (key) and URI (value) of all necessary XML namespaces.
    scene_id: str
        Scene basename.
    extent: str
        Footprint coordinates of the scene.
    center: str
        Center coordinates of the footprint.
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
