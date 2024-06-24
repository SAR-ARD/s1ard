import os
import re
import sys
import shutil
from statistics import mean
from datetime import datetime, timezone
import pystac
from pystac.extensions.sar import SarExtension, FrequencyBand, Polarization, ObservationDirection
from pystac.extensions.sat import SatExtension, OrbitState
from pystac.extensions.projection import ProjectionExtension
from pystac.extensions.view import ViewExtension
from pystac.extensions.mgrs import MgrsExtension
from pystac.extensions.file import FileExtension, ByteOrder
from pystac.extensions.raster import RasterExtension, RasterBand, DataType
from pystac.extensions.classification import ClassificationExtension, Classification
from spatialist import Raster
from spatialist.ancillary import finder
from s1ard.metadata.mapping import ASSET_MAP
from s1ard.metadata.extract import get_header_size
from s1ard.ancillary import compute_hash
import logging

log = logging.getLogger('s1ard')


def parse(meta, target, assets, exist_ok=False):
    """
    Wrapper for :func:`~s1ard.metadata.stac.source_json` and :func:`~s1ard.metadata.stac.product_json`.
    
    Parameters
    ----------
    meta: dict
        Metadata dictionary generated with :func:`~s1ard.metadata.extract.meta_dict`
    target: str
        A path pointing to the root directory of a product scene.
    assets: list[str]
        List of paths to all GeoTIFF and VRT assets of the currently processed ARD product.
    exist_ok: bool
        Do not create files if they already exist?
    """
    source_json(meta=meta, target=target, exist_ok=exist_ok)
    product_json(meta=meta, target=target, assets=assets, exist_ok=exist_ok)


def source_json(meta, target, exist_ok=False):
    """
    Function to generate source-level metadata for an ARD product in STAC compliant JSON format.
    
    Parameters
    ----------
    meta: dict
        Metadata dictionary generated with :func:`~s1ard.metadata.extract.meta_dict`.
    target: str
        A path pointing to the root directory of a product scene.
    exist_ok: bool
        Do not create files if they already exist?
    """
    metadir = os.path.join(target, 'source')
    for uid in list(meta['source'].keys()):
        scene = os.path.basename(meta['source'][uid]['filename']).split('.')[0]
        outname = os.path.join(metadir, '{}.json'.format(scene))
        if os.path.isfile(outname) and exist_ok:
            continue
        log.info(f'creating {outname}')
        start = meta['source'][uid]['timeStart']
        stop = meta['source'][uid]['timeStop']
        date = start + (stop - start) / 2
        
        # Initialise STAC item
        item = pystac.Item(id=scene,
                           geometry=meta['source'][uid]['geom_stac_geometry_4326'],
                           bbox=meta['source'][uid]['geom_stac_bbox_4326'],
                           datetime=date,
                           properties={})
        
        # Add common metadata
        item.common_metadata.start_datetime = start
        item.common_metadata.end_datetime = stop
        item.common_metadata.created = datetime.strptime(meta['source'][uid]['processingDate'], '%Y-%m-%dT%H:%M:%S.%f')
        item.common_metadata.instruments = [meta['common']['instrumentShortName'].lower()]
        item.common_metadata.constellation = meta['common']['constellation']
        item.common_metadata.platform = meta['common']['platformFullname']
        
        # Initialise STAC extensions for properties
        sar_ext = SarExtension.ext(item, add_if_missing=True)
        sat_ext = SatExtension.ext(item, add_if_missing=True)
        view_ext = ViewExtension.ext(item, add_if_missing=True)
        item.stac_extensions.append('https://stac-extensions.github.io/processing/v1.1.0/schema.json')
        item.stac_extensions.append('https://stac-extensions.github.io/card4l/v0.1.0/sar/source.json')
        
        # Add properties
        sat_ext.apply(orbit_state=OrbitState[meta['common']['orbitDirection'].upper()],
                      relative_orbit=meta['common']['orbitNumber_rel'],
                      absolute_orbit=meta['common']['orbitNumber_abs'],
                      anx_datetime=datetime.strptime(meta['source'][uid]['ascendingNodeDate'], '%Y-%m-%dT%H:%M:%S.%f'))
        sar_ext.apply(instrument_mode=meta['common']['operationalMode'],
                      frequency_band=FrequencyBand[meta['common']['radarBand'].upper()],
                      polarizations=[Polarization[pol] for pol in meta['common']['polarisationChannels']],
                      product_type=meta['source'][uid]['productType'],
                      center_frequency=float(meta['common']['radarCenterFreq'] / 1e9),
                      resolution_range=min(meta['source'][uid]['rangeResolution'].values()),
                      resolution_azimuth=min(meta['source'][uid]['azimuthResolution'].values()),
                      pixel_spacing_range=mean(meta['source'][uid]['rangePixelSpacing'].values()),
                      pixel_spacing_azimuth=mean(meta['source'][uid]['azimuthPixelSpacing'].values()),
                      looks_range=mean(meta['source'][uid]['rangeNumberOfLooks'].values()),
                      looks_azimuth=mean(meta['source'][uid]['azimuthNumberOfLooks'].values()),
                      looks_equivalent_number=meta['source'][uid]['perfEquivalentNumberOfLooks'],
                      observation_direction=ObservationDirection[meta['common']['antennaLookDirection']])
        view_ext.apply(incidence_angle=meta['source'][uid]['incidenceAngleMidSwath'],
                       azimuth=meta['source'][uid]['instrumentAzimuthAngle'])
        item.properties['processing:facility'] = meta['source'][uid]['processingCenter']
        item.properties['processing:software'] = {meta['source'][uid]['processorName']:
                                                      meta['source'][uid]['processorVersion']}
        item.properties['processing:level'] = meta['common']['processingLevel']
        item.properties['card4l:specification'] = meta['prod']['productName-short']
        item.properties['card4l:specification_version'] = meta['prod']['card4l-version']
        item.properties['card4l:beam_id'] = meta['common']['swathIdentifier']
        item.properties['card4l:orbit_data_source'] = meta['source'][uid]['orbitDataSource']
        item.properties['card4l:orbit_mean_altitude'] = float(meta['common']['orbitMeanAltitude'])
        range_look_bandwidth = {k: v / 1e9 for k, v in meta['source'][uid]['rangeLookBandwidth'].items()}  # GHz
        azimuth_look_bandwidth = {k: v / 1e9 for k, v in meta['source'][uid]['azimuthLookBandwidth'].items()}  # GHz
        item.properties['card4l:source_processing_parameters'] = {'lut_applied': meta['source'][uid]['lutApplied'],
                                                                  'range_look_bandwidth': range_look_bandwidth,
                                                                  'azimuth_look_bandwidth': azimuth_look_bandwidth}
        for field, key in zip(['card4l:resolution_range', 'card4l:resolution_azimuth'],
                              ['rangeResolution', 'azimuthResolution']):
            res = {}
            for k, v in meta['source'][uid][key].items():
                res[k] = float(v)
            item.properties[field] = res
        item.properties['card4l:source_geometry'] = meta['source'][uid]['dataGeometry']
        item.properties['card4l:incidence_angle_near_range'] = meta['source'][uid]['incidenceAngleMin']
        item.properties['card4l:incidence_angle_far_range'] = meta['source'][uid]['incidenceAngleMax']
        item.properties['card4l:noise_equivalent_intensity'] = meta['source'][uid]['perfEstimates']
        item.properties['card4l:noise_equivalent_intensity_type'] = meta['source'][uid][
            'perfNoiseEquivalentIntensityType']
        item.properties['card4l:peak_sidelobe_ratio'] = meta['source'][uid]['perfPeakSideLobeRatio']
        item.properties['card4l:integrated_sidelobe_ratio'] = meta['source'][uid]['perfIntegratedSideLobeRatio']
        item.properties['card4l:mean_faraday_rotation_angle'] = meta['source'][uid]['faradayMeanRotationAngle']
        item.properties['card4l:ionosphere_indicator'] = meta['source'][uid]['ionosphereIndicator']
        
        # Add links
        item.add_link(link=pystac.Link(rel='card4l-document',
                                       target=meta['prod']['card4l-link'].replace('.pdf', '.docx'),
                                       media_type='application/vnd.openxmlformats-officedocument.wordprocessingml'
                                                  '.document',
                                       title='CARD4L Product Family Specification: {} (v{})'
                                             ''.format(meta['prod']['productName'], meta['prod']['card4l-version'])))
        item.add_link(link=pystac.Link(rel='card4l-document',
                                       target=meta['prod']['card4l-link'],
                                       media_type='application/pdf',
                                       title='CARD4L Product Family Specification: {} (v{})'
                                             ''.format(meta['prod']['productName'], meta['prod']['card4l-version'])))
        item.add_link(link=pystac.Link(rel='about',
                                       target=meta['source'][uid]['doi'],
                                       title='Product definition reference.'))
        item.add_link(link=pystac.Link(rel='access',
                                       target=meta['source'][uid]['access'],
                                       title='Product data access.'))
        item.add_link(link=pystac.Link(rel='satellite',
                                       target=meta['common']['platformReference'],
                                       title='CEOS Missions, Instruments and Measurements Database record'))
        item.add_link(link=pystac.Link(rel='state-vectors',
                                       target=meta['source'][uid]['orbitStateVector'],
                                       title='Orbit data file containing state vectors.'))
        item.add_link(link=pystac.Link(rel='sensor-calibration',
                                       target=meta['source'][uid]['sensorCalibration'],
                                       title='Reference describing sensor calibration parameters.'))
        item.add_link(link=pystac.Link(rel='pol-cal-matrices',
                                       target=meta['source'][uid]['polCalMatrices'],
                                       title='Reference to the complex-valued polarimetric distortion matrices.'))
        item.add_link(link=pystac.Link(rel='referenced-faraday-rotation',
                                       target=meta['source'][uid]['faradayRotationReference'],
                                       title='Reference describing the method used to derive the estimate for the mean'
                                             ' Faraday rotation angle.'))
        
        # Add assets
        xml_relpath = './' + os.path.relpath(outname.replace('.json', '.xml'), metadir).replace('\\', '/')
        item.add_asset(key='card4l',
                       asset=pystac.Asset(href=xml_relpath,
                                          title='Metadata in XML format.',
                                          media_type=pystac.MediaType.XML,
                                          roles=['metadata', 'card4l']))
        _asset_add_orig_src(metadir=metadir, uid=uid, item=item)
        
        item.save_object(dest_href=outname)


def _asset_add_orig_src(metadir, uid, item):
    """
    Helper function to add the original source metadata files as assets to a STAC item.
    
    Parameters
    ----------
    metadir: str
        Source directory of the current ARD product.
    uid: str
        Unique identifier of a source scene.
    item: pystac.Item
        The pystac.Item to add the assets to.
    
    Returns
    -------
    
    """
    pattern = r'^(.+?)-\d{8}t\d{6}-\d{8}t\d{6}-\w+-\w+-\d{3}\.xml$'
    prefixes = {'calibration': 'Calibration metadata',
                'noise': 'Estimated thermal noise look-up tables',
                'rfi': 'Radio Frequency Interference metadata'}
    
    root_dir = os.path.join(metadir, uid)
    if not os.path.isdir(root_dir):
        return
    file_list = finder(target=root_dir, matchlist=['*.safe', '*.xml'], foldermode=0)
    if len(file_list) > 0:
        for file in file_list:
            basename = os.path.basename(file)
            href = './' + os.path.relpath(file, metadir).replace('\\', '/')
            if basename == 'manifest.safe':
                key = 'manifest'
                title = 'Mandatory product metadata'
            else:
                try:
                    key = re.match(pattern, basename).group(1)
                except AttributeError:
                    raise RuntimeError(
                        'Unexpected file in original source metadata directory: ' + os.path.join(root_dir, file))
                title = prefixes.get(key.split('-')[0], 'Measurement metadata')
                if title == 'Measurement metadata':
                    title = title + ' ({},{})'.format(key.split('-')[1].upper(), key.split('-')[3].upper())
                else:
                    title = title + ' ({},{})'.format(key.split('-')[2].upper(), key.split('-')[4].upper())
            
            item.add_asset(key=key,
                           asset=pystac.Asset(href=href,
                                              title=title,
                                              media_type=pystac.MediaType.XML,
                                              roles=['metadata']))


def product_json(meta, target, assets, exist_ok=False):
    """
    Function to generate product-level metadata for an ARD product in STAC compliant JSON format.
    
    Parameters
    ----------
    meta: dict
        Metadata dictionary generated with :func:`~s1ard.metadata.extract.meta_dict`.
    target: str
        A path pointing to the root directory of a product scene.
    assets: list[str]
        List of paths to all GeoTIFF and VRT assets of the currently processed ARD product.
    exist_ok: bool
        Do not create files if they already exist?
    """
    scene_id = os.path.basename(target)
    outname = os.path.join(target, '{}.json'.format(scene_id))
    if os.path.isfile(outname) and exist_ok:
        return
    log.info(f'creating {outname}')
    start = meta['prod']['timeStart']
    stop = meta['prod']['timeStop']
    date = start + (stop - start) / 2
    mgrs = meta['prod']['mgrsID']
    
    # Initialise STAC item
    item = pystac.Item(id=scene_id,
                       geometry=meta['prod']['geom_stac_geometry_4326'],
                       bbox=meta['prod']['geom_stac_bbox_4326'],
                       datetime=date,
                       properties={})
    
    # Add common metadata
    item.common_metadata.license = meta['prod']['licence']
    item.common_metadata.start_datetime = start
    item.common_metadata.end_datetime = stop
    item.common_metadata.created = meta['prod']['timeCreated']
    item.common_metadata.instruments = [meta['common']['instrumentShortName'].lower()]
    item.common_metadata.constellation = meta['common']['constellation']
    item.common_metadata.platform = meta['common']['platformFullname']
    item.common_metadata.gsd = float(meta['prod']['pxSpacingColumn'])
    
    # Initialise STAC extensions for properties
    sar_ext = SarExtension.ext(item, add_if_missing=True)
    sat_ext = SatExtension.ext(item, add_if_missing=True)
    proj_ext = ProjectionExtension.ext(item, add_if_missing=True)
    mgrs_ext = MgrsExtension.ext(item, add_if_missing=True)
    item.stac_extensions.append('https://stac-extensions.github.io/processing/v1.1.0/schema.json')
    item.stac_extensions.append('https://stac-extensions.github.io/card4l/v0.1.0/sar/product.json')
    
    # Add properties
    sat_ext.apply(orbit_state=OrbitState[meta['common']['orbitDirection'].upper()],
                  relative_orbit=meta['common']['orbitNumber_rel'],
                  absolute_orbit=meta['common']['orbitNumber_abs'])
    sar_ext.apply(instrument_mode=meta['common']['operationalMode'],
                  frequency_band=FrequencyBand[meta['common']['radarBand'].upper()],
                  polarizations=[Polarization[pol] for pol in meta['common']['polarisationChannels']],
                  product_type=meta['prod']['productName-short'],
                  looks_range=meta['prod']['rangeNumberOfLooks'],
                  looks_azimuth=meta['prod']['azimuthNumberOfLooks'],
                  looks_equivalent_number=meta['prod']['equivalentNumberLooks'])
    proj_ext.apply(epsg=int(meta['prod']['crsEPSG']),
                   wkt2=meta['prod']['crsWKT'],
                   bbox=meta['prod']['geom_stac_bbox_native'],
                   shape=[int(meta['prod']['numPixelsPerLine']), int(meta['prod']['numLines'])],
                   transform=meta['prod']['transform'])
    mgrs_ext.apply(latitude_band=mgrs[2:3],
                   grid_square=mgrs[3:],
                   utm_zone=int(mgrs[:2]))
    if meta['prod']['processingCenter'] is not None:
        item.properties['processing:facility'] = meta['prod']['processingCenter']
    item.properties['processing:software'] = {meta['prod']['processorName']: meta['prod']['processorVersion']}
    item.properties['processing:level'] = meta['common']['processingLevel']
    item.properties['card4l:specification'] = meta['prod']['productName-short']
    item.properties['card4l:specification_version'] = meta['prod']['card4l-version']
    item.properties['card4l:beam_id'] = meta['common']['swathIdentifier']
    item.properties['card4l:measurement_type'] = meta['prod']['backscatterMeasurement']
    item.properties['card4l:measurement_convention'] = meta['prod']['backscatterConvention']
    item.properties['card4l:pixel_coordinate_convention'] = meta['prod']['pixelCoordinateConvention']
    item.properties['card4l:speckle_filtering'] = meta['prod']['speckleFilterApplied']
    item.properties['card4l:noise_removal_applied'] = meta['prod']['noiseRemovalApplied']
    item.properties['card4l:conversion_eq'] = meta['prod']['backscatterConversionEq']
    item.properties['card4l:relative_radiometric_accuracy'] = meta['prod']['radiometricAccuracyRelative']
    item.properties['card4l:absolute_radiometric_accuracy'] = meta['prod']['radiometricAccuracyAbsolute']
    item.properties['card4l:resampling_method'] = meta['prod']['geoCorrResamplingMethod']
    item.properties['card4l:dem_resampling_method'] = meta['prod']['demResamplingMethod']
    item.properties['card4l:egm_resampling_method'] = meta['prod']['demEGMResamplingMethod']
    item.properties['card4l:geometric_accuracy_type'] = meta['prod']['geoCorrAccuracyType']
    for x in ['Northern', 'Eastern']:
        key = ['geoCorrAccuracy{}{}'.format(x, y) for y in ['STDev', 'Bias']]
        stddev = float(meta['prod'][key[0]]) if meta['prod'][key[0]] is not None else None
        bias = float(meta['prod'][key[1]]) if meta['prod'][key[1]] is not None else None
        item.properties['card4l:{}_geometric_accuracy'.format(x.lower())] = {'bias': bias, 'stddev': stddev}
    item.properties['card4l:geometric_accuracy_radial_rmse'] = meta['prod']['geoCorrAccuracy_rRMSE']
    
    # Add links
    item.add_link(link=pystac.Link(rel='card4l-document',
                                   target=meta['prod']['card4l-link'].replace('.pdf', '.docx'),
                                   media_type='application/vnd.openxmlformats-officedocument.wordprocessingml.document',
                                   title='CARD4L Product Family Specification: {} (v{})'
                                         ''.format(meta['prod']['productName'], meta['prod']['card4l-version'])))
    item.add_link(link=pystac.Link(rel='card4l-document',
                                   target=meta['prod']['card4l-link'],
                                   media_type='application/pdf',
                                   title='CARD4L Product Family Specification: {} (v{})'
                                         ''.format(meta['prod']['productName'], meta['prod']['card4l-version'])))
    for src in list(meta['source'].keys()):
        x = os.path.basename(meta['source'][src]['filename']).split('.')[0]
        src_target = os.path.join('./source', '{}.json'.format(x)).replace('\\', '/')
        item.add_link(link=pystac.Link(rel='derived_from',
                                       target=src_target,
                                       media_type='application/json',
                                       title='Source metadata formatted in STAC compliant JSON format.'))
    if meta['prod']['doi'] is not None:
        item.add_link(link=pystac.Link(rel='about',
                                       target=meta['prod']['doi'],
                                       title='Product definition reference.'))
    if meta['prod']['access'] is not None:
        item.add_link(link=pystac.Link(rel='access',
                                       target=meta['prod']['access'],
                                       title='Product data access.'))
    item.add_link(link=pystac.Link(rel='related',
                                   target=meta['prod']['ancillaryData_KML'],
                                   title='Sentinel-2 Military Grid Reference System (MGRS) tiling grid file '
                                         'used as auxiliary data during processing.'))
    if meta['prod']['noiseRemovalApplied']:
        item.add_link(link=pystac.Link(rel='noise-removal',
                                       target=meta['prod']['noiseRemovalAlgorithm'],
                                       title='Reference to the noise removal algorithm details.'))
    if meta['prod']['RTCAlgorithm'] is not None:
        item.add_link(link=pystac.Link(rel='radiometric-terrain-correction',
                                       target=meta['prod']['RTCAlgorithm'],
                                       title='Reference to the Radiometric Terrain Correction algorithm details.'))
    if meta['prod']['windNormReferenceModel'] is not None:
        item.add_link(link=pystac.Link(rel='wind-norm-reference',
                                       target=meta['prod']['windNormReferenceModel'],
                                       title='Reference to the model used to create the wind normalisation layer.'))
    item.add_link(link=pystac.Link(rel='radiometric-accuracy',
                                   target=meta['prod']['radiometricAccuracyReference'],
                                   title='Reference describing the radiometric uncertainty of the product.'))
    item.add_link(link=pystac.Link(rel='geometric-correction',
                                   target=meta['prod']['geoCorrAlgorithm'],
                                   title='Reference to the Geometric Correction algorithm details.'))
    item.add_link(link=pystac.Link(rel='{}-model'.format(meta['prod']['demType']),
                                   target=meta['prod']['demReference'],
                                   title='Digital Elevation Model used as auxiliary data during processing: '
                                         '{}'.format(meta['prod']['demName'])))
    item.add_link(link=pystac.Link(rel='earth-gravitational-model',
                                   target=meta['prod']['demEGMReference'],
                                   title='Reference to the Earth Gravitational Model (EGM) used for Geometric '
                                         'Correction.'))
    item.add_link(link=pystac.Link(rel='geometric-accuracy',
                                   target=meta['prod']['geoCorrAccuracyReference'],
                                   title='Reference documenting the estimate of absolute localization error.'))
    item.add_link(link=pystac.Link(rel='gridding-convention',
                                   target=meta['prod']['griddingConventionURL'],
                                   title='Reference describing the gridding convention used.'))
    
    # Add assets
    assets = assets.copy()
    assets.append(outname.replace('.json', '.xml'))
    
    assets_dict = {'measurement': {},
                   'annotation': {},
                   'metadata': {}}
    for asset in assets:
        relpath = './' + os.path.relpath(asset, target).replace('\\', '/')
        
        size = os.path.getsize(asset)
        hash = compute_hash(asset)
        created = None
        header_size = None
        media_type = pystac.MediaType.XML  # VRT
        byte_order = None
        if asset.endswith('.tif'):
            with Raster(asset) as ras:
                nodata = ras.nodata
            created = datetime.fromtimestamp(os.path.getctime(asset), tz=timezone.utc).isoformat()
            header_size = get_header_size(tif=asset)
            media_type = pystac.MediaType.COG
            byte_order = ByteOrder.LITTLE_ENDIAN
        
        if 'measurement' in asset:
            key, title = _asset_get_key_title(meta=meta, asset=asset)
            stac_asset = pystac.Asset(href=relpath,
                                      title=title,
                                      media_type=media_type,
                                      roles=['backscatter', 'data'],
                                      extra_fields=None)
            if asset.endswith('.tif'):
                stac_asset.extra_fields = {'created': created,
                                           'card4l:border_pixels': meta['prod']['numBorderPixels']}
                _asset_handle_raster_ext(stac_asset=stac_asset, nodata=nodata)
            file_ext = FileExtension.ext(stac_asset)
            file_ext.apply(byte_order=byte_order, size=size,
                           header_size=header_size, checksum=hash)
            assets_dict['measurement'][key] = stac_asset
        elif 'annotation' in asset:
            key, title = _asset_get_key_title(meta=meta, asset=asset)
            if key == '-np-[vh]{2}.tif':
                pol = re.search('-[vh]{2}', relpath).group().removeprefix('-')
                asset_key = 'noise-power-{}'.format(pol)
                title = "{} {}".format(title, pol.upper())
            else:
                asset_key = ASSET_MAP[key]['role']
            
            if ASSET_MAP[key]['unit'] is None:
                ASSET_MAP[key]['unit'] = 'unitless'
            
            stac_asset = pystac.Asset(href=relpath,
                                      title=title,
                                      media_type=media_type,
                                      roles=[ASSET_MAP[key]['role'], 'metadata'],
                                      extra_fields=None)
            if key == '-ei.tif':
                stac_asset.extra_fields = {'card4l:ellipsoidal_height': meta['prod']['ellipsoidalHeight']}
            file_ext = FileExtension.ext(stac_asset)
            file_ext.apply(byte_order=byte_order, size=size,
                           header_size=header_size, checksum=hash)
            _asset_handle_raster_ext(stac_asset=stac_asset, nodata=nodata, key=key, meta=meta, asset=asset)
            assets_dict['annotation'][asset_key] = stac_asset
        else:
            stac_asset = pystac.Asset(href=relpath,
                                      title='Metadata in XML format.',
                                      media_type=pystac.MediaType.XML,
                                      roles=['metadata', 'card4l'])
            file_ext = FileExtension.ext(stac_asset)
            file_ext.apply(byte_order=byte_order, size=size,
                           header_size=header_size, checksum=hash)
            assets_dict['metadata']['card4l'] = stac_asset
    
    for category in ['measurement', 'annotation', 'metadata']:
        for key in sorted(assets_dict[category].keys()):
            item.add_asset(key=key, asset=assets_dict[category][key])
    
    if any(x in item.get_assets().keys() for x in ['acquisition-id', 'data-mask']):
        ClassificationExtension.add_to(item)
    FileExtension.add_to(item)
    RasterExtension.add_to(item)
    item.save_object(dest_href=outname)


def _asset_get_key_title(meta, asset):
    """
    Helper function to handle creation of identifying key and title of a given asset.
    
    Parameters
    ----------
    meta: dict
        Metadata dictionary generated with :func:`~s1ard.metadata.extract.meta_dict`.
    asset: str
        Path to a GeoTIFF or VRT asset.
    
    Returns
    -------
    key: str
        Key identifying the asset.
    title: str
        Generated title of the asset.
    """
    key = None
    title = None
    if 'measurement' in asset:
        title_dict = {'g': 'gamma nought', 's': 'sigma nought',
                      'lin': 'linear', 'log': 'logarithmic'}
        pattern = ('(?P<key>(?P<pol>[vhc]{2})-'
                   '(?P<nought>[gs])-'
                   '(?P<scaling>lin|log)'
                   '(?P<windnorm>-wn|))')
        info = re.search(pattern, asset).groupdict()
        key = info['key']
        
        if re.search('cc-[gs]-lin', key):
            pols = meta['common']['polarisationChannels']
            co = pols.pop(0) if pols[0][0] == pols[0][1] else pols.pop(1)
            cross = pols[0]
            title = 'RGB color composite (' \
                    '{co}-{nought}-lin, ' \
                    '{cross}-{nought}-lin, ' \
                    '{co}-{nought}-lin/{cross}-{nought}-lin)'
            title = title.format(co=co.lower(),
                                 cross=cross.lower(),
                                 nought=info['nought'])
        else:
            if info['windnorm'] == '-wn':
                skeleton = '{pol} {nought} {subtype} wind normalisation ratio, {scale} scaling'
            else:
                skeleton = '{pol} {nought} {subtype} backscatter, {scale} scaling'
            subtype = 'RTC'
            title = skeleton.format(pol=info['pol'].upper(),
                                    nought=title_dict[info['nought']],
                                    subtype=subtype,
                                    scale=title_dict[info['scaling']])
    elif 'annotation' in asset:
        key = re.search('-[a-z]{2}(?:-[a-z]{2}|).tif', asset).group()
        np_pat = '-np-[vh]{2}.tif'
        if re.search(np_pat, key) is not None:
            key = np_pat
        title = ASSET_MAP[key]['title']
    return key, title


def _asset_handle_raster_ext(stac_asset, nodata, key=None, meta=None, asset=None):
    """
    Helper function to handle the STAC RasterExtension for a given pystac.Asset
    
    Parameters
    ----------
    stac_asset: pystac.Asset
        The pystac.Asset to set the RasterExtension for.
    nodata: float
        The asset's nodata value.
    key: str
        Key identifying the asset. Only necessary for annotation assets.
    meta: dict
        Metadata dictionary generated with :func:`~s1ard.metadata.extract.meta_dict`.
        Only necessary for annotation assets.
    asset: str
        Path to a GeoTIFF or VRT asset. Only necessary for annotation assets.
    
    Returns
    -------
    None
    """
    raster_ext = RasterExtension.ext(stac_asset)
    if 'measurement' in stac_asset.href:
        raster_ext.bands = [RasterBand.create(nodata=nodata,
                                              data_type=DataType.FLOAT32,
                                              unit='natural')]
    elif 'annotation' in stac_asset.href:
        if key is None or meta is None or asset is None:
            raise ValueError(f'`key`, `meta` and `asset` parameters need to be defined to handle RasterExtension '
                             f'for {stac_asset.href}')
        if key == '-id.tif':
            src_list = [
                os.path.basename(meta['source'][src]['filename']).replace('.SAFE', '').replace('.zip', '')
                for src in list(meta['source'].keys())]
            band = RasterBand.create(nodata=nodata,
                                     data_type=DataType.UINT8,
                                     unit=ASSET_MAP[key]['unit'])
            class_ext = ClassificationExtension.ext(band)
            class_ext.classes = [Classification.create(value=i + 1, description=j) for i, j in enumerate(src_list)]
            raster_ext.bands = [band]
        elif key == '-dm.tif':
            with Raster(asset) as dm_ras:
                band_descr = [dm_ras.raster.GetRasterBand(band).GetDescription() for band in
                              range(1, dm_ras.bands + 1)]
            samples = [x for x in band_descr if x in ASSET_MAP[key]['allowed']]
            bands = []
            for sample in samples:
                band = RasterBand.create(nodata=nodata,
                                         data_type=DataType.UINT8,
                                         unit=ASSET_MAP[key]['unit'])
                class_ext = ClassificationExtension.ext(band, add_if_missing=True)
                class_ext.classes = [Classification.create(value=1, description=sample)]
                bands.append(band)
            raster_ext.bands = bands
        else:
            raster_ext.bands = [RasterBand.create(nodata=nodata,
                                                  data_type=DataType.FLOAT32,
                                                  unit=ASSET_MAP[key]['unit'])]
    if key == '-em.tif':
        raster_ext.bands[0].spatial_resolution = int(meta['prod']['demGSD'].split()[0])


def make_catalog(directory, product_type, recursive=True, silent=False):
    """
    For a given directory of Sentinel-1 ARD products, this function will create a high-level STAC
    :class:`~pystac.catalog.Catalog` object serving as the STAC endpoint and lower-level STAC
    :class:`~pystac.collection.Collection` objects for each subdirectory corresponding to a unique MGRS tile ID.
    
    WARNING: The directory content will be reorganized into subdirectories based on the ARD type and unique MGRS tile
    IDs if this is not yet the case.
    
    Parameters
    ----------
    directory: str
        Path to a directory that contains ARD products.
    product_type: str
        Type of ARD products. Options: 'NRB' or 'ORB'.
    recursive: bool, optional
        Search `directory` recursively? Default is True.
    silent: bool, optional
        Should the output during directory reorganization be suppressed? Default is False.
    
    Returns
    -------
    nrb_catalog: pystac.catalog.Catalog
        STAC Catalog object
    
    Notes
    -----
    The returned STAC Catalog object contains Item asset hrefs that are absolute, whereas the actual on-disk files
    contain relative asset hrefs corresponding to the self-contained Catalog-Type. The returned in-memory STAC Catalog
    object deviates in this regard to ensure compatibility with the stackstac library:
    https://github.com/gjoseph92/stackstac/issues/20
    """
    overwrite = False
    product_type = product_type.upper()
    pattern = fr'^S1[AB]_(IW|EW|S[1-6])_{product_type}__1S(SH|SV|DH|DV|VV|HH|HV|VH)_[0-9]{{8}}T[0-9]{{6}}_[0-9]{{6}}_' \
              fr'[0-9A-F]{{6}}_[0-9A-Z]{{5}}_[0-9A-Z]{{4}}$'
    products = finder(target=directory, matchlist=[pattern], foldermode=2, regex=True, recursive=recursive)
    directory = os.path.join(directory, product_type)
    
    # Check if Catalog already exists
    catalog_path = os.path.join(directory, 'catalog.json')
    if os.path.isfile(catalog_path):
        overwrite = True
        catalog = pystac.Catalog.from_file(catalog_path)
        items = catalog.get_all_items()
        item_ids = [item.id for item in items]
        products_base = [os.path.basename(prod) for prod in products]
        diff = set(products_base) - set(item_ids)
        if len(diff) == 0:
            # See note in docstring - https://github.com/gjoseph92/stackstac/issues/20
            catalog.make_all_asset_hrefs_absolute()
            log.info(f"existing STAC endpoint found: {os.path.join(directory, 'catalog.json')}")
            return catalog
    
    sp_extent = pystac.SpatialExtent([None, None, None, None])
    tmp_extent = pystac.TemporalExtent([None, None])
    
    unique_tiles = list(
        set([re.search(re.compile(r'_[0-9A-Z]{5}_'), prod).group().replace('_', '') for prod in products]))
    products = _reorganize_by_tile(directory=directory, product_type=product_type, products=products,
                                   recursive=recursive, silent=silent)
    
    catalog = pystac.Catalog(id=f'{product_type.lower()}_catalog',
                             description=f'STAC Catalog of Sentinel-1 {product_type} products.',
                             title=f'STAC Catalog of Sentinel-1 {product_type} products.',
                             catalog_type=pystac.CatalogType.SELF_CONTAINED)
    
    for tile in unique_tiles:
        tile_collection = pystac.Collection(id=tile,
                                            description=f'STAC Collection of Sentinel-1 {product_type} products for '
                                                        f'MGRS tile {tile}.',
                                            title=f'STAC Collection of Sentinel-1 {product_type} products for '
                                                  f'MGRS tile {tile}.',
                                            extent=pystac.Extent(sp_extent, tmp_extent),
                                            keywords=['sar', 'backscatter', 'esa', 'copernicus', 'sentinel'],
                                            providers=[pystac.Provider(name='ESA',
                                                                       roles=[pystac.ProviderRole.LICENSOR,
                                                                              pystac.ProviderRole.PRODUCER])])
        catalog.add_child(tile_collection)
        
        items = []
        for prod in products:
            if tile in prod:
                item_path = os.path.join(prod, os.path.basename(prod) + '.json')
                item = pystac.read_file(href=item_path)
                items.append(item)
                tile_collection.add_item(item=item)
            else:
                continue
        
        extent = tile_collection.extent.from_items(items=items)
        tile_collection.extent = extent
    
    # Save Catalog and Collections on disk
    catalog.normalize_and_save(root_href=directory)
    
    # See note in docstring - https://github.com/gjoseph92/stackstac/issues/20
    catalog.make_all_asset_hrefs_absolute()
    
    if overwrite:
        log.info(f"existing STAC endpoint updated: {os.path.join(directory, 'catalog.json')}")
    else:
        log.info(f"new STAC endpoint created: {os.path.join(directory, 'catalog.json')}")
    return catalog


def _reorganize_by_tile(directory, product_type, products=None, recursive=True, silent=False):
    """
    Reorganizes a directory containing Sentinel-1 ARD products based on the ARD type and unique MGRS tile IDs.
    
    Parameters
    ----------
    directory: str
        Path to a directory that contains ARD products.
    product_type: str
        Type of ARD products. Options: 'NRB' or 'ORB'.
    products: list[str] or None, optional
        List of ARD product paths. Will be created from `directory` if not provided.
    recursive: bool, optional
        Search `directory` recursively? Default is True.
    silent: bool, optional
        If False (default), a message for each ARD product is printed if it has been moved to a new location or not.
    
    Returns
    -------
    products_new: list[str]
        An updated list of ARD product paths.
    """
    if products is None:
        parent_dir = os.path.dirname(directory)
        pattern = fr'^S1[AB]_(IW|EW|S[1-6])_{product_type}__1S(SH|SV|DH|DV|VV|HH|HV|VH)_[0-9]{{8}}T[0-9]{{6}}_' \
                  fr'[0-9]{{6}}_[0-9A-F]{{6}}_[0-9A-Z]{{5}}_[0-9A-Z]{{4}}$'
        products = finder(target=parent_dir, matchlist=[pattern], foldermode=2, regex=True, recursive=recursive)
    
    inp = input('WARNING:\n{}\nand the ARD products it contains will be reorganized into subdirectories '
                'based on unique MGRS tile IDs if this directory structure does not yet exist. '
                '\nDo you wish to continue? [yes|no] '.format(directory))
    if inp == 'yes':
        tile_dict = {}
        for prod in products:
            tile = re.search(re.compile(r'_[0-9A-Z]{5}_'), prod).group().replace('_', '')
            if tile in tile_dict and isinstance(tile_dict[tile], list):
                tile_dict[tile].append(prod)
            else:
                tile_dict[tile] = [prod]
        
        tiles = list(tile_dict.keys())
        products_new = []
        for tile in tiles:
            tile_dir = os.path.join(directory, tile)
            os.makedirs(tile_dir, exist_ok=True)
            
            for old_dir in tile_dict[tile]:
                new_dir = os.path.join(tile_dir, os.path.basename(old_dir))
                products_new.append(new_dir)
                
                if os.path.dirname(old_dir) != tile_dir:
                    shutil.move(old_dir, new_dir)
                    if not silent:
                        log.info(f"-> {os.path.basename(old_dir)} moved to {tile_dir}")
                else:
                    if not silent:
                        log.info(f"xx {os.path.basename(old_dir)} already in {tile_dir} (skip!)")
                    continue
        return products_new
    else:
        log.info('abort!')
        sys.exit(0)
