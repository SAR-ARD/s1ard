import os
import re
import shutil
import zipfile
import math
from typing import Union, Literal
from collections import defaultdict
from statistics import mean
from lxml import etree
from dateutil.parser import parse as dateparse
from datetime import timezone, datetime
import numpy as np
from spatialist import Raster
from spatialist.auxil import gdalwarp
from spatialist.ancillary import finder, dissolve
from spatialist.raster import rasterize
from pyroSAR.drivers import ID
from pyroSAR.ancillary import get_geometry
from osgeo import gdal
from s1ard.config import version_dict
from s1ard.metadata.mapping import (RES_MAP_SLC, RES_MAP_GRD, ENL_MAP_GRD,
                                    OSV_MAP, SLC_ACC_MAP, URL)
from s1ard.processors.registry import load_processor
from cesard.ancillary import get_tmp_name, defaultdict_to_dict
from cesard.metadata.extract import (calc_enl, calc_performance_estimates,
                                     geometry_from_vec)
from cesard.metadata.mapping import DEM_MAP, LERC_ERR_THRES
from cesard.metadata.model import (
    ARDMetadata, AxisAccuracy, BackscatterMetadata, Card4lMetadata,
    CommonMetadata, CompressionMetadata, DEMMetadata, GeometricAccuracyMetadata,
    GeometricCorrectionMetadata, GeometryMetadata, GridMetadata,
    GroundSamplingDistance, IncidenceAngleMetadata, NoiseRemovalMetadata,
    NOT_IMPLEMENTED_NUMBER, NOT_IMPLEMENTED_TEXT, PerformanceEstimate,
    ProcessingMetadata, ProductMetadata, RadiometricAccuracyMetadata,
    SourceMetadata, SourceOrbitMetadata, SourcePerformanceMetadata,
    SwathAxisMetadata, WindNormalizationMetadata,
)
from typing import Any

gdal.UseExceptions()


def append_wind_norm(
        meta: ARDMetadata,
        wm_ref_speed: list[str] | None,
        wm_ref_direction: list[str] | None
) -> None:
    """
    Add wind model information to an ARD metadata object.

    Parameters
    ----------
    meta
        Metadata extracted by :func:`meta_dict`.
    wm_ref_speed
        List of paths pointing to the wind model reference speed files.
    wm_ref_direction
        List of paths pointing to the wind model reference direction files.
    """
    if wm_ref_speed is not None and wm_ref_direction is not None:
        wm_ref_mean_speed, wm_ref_mean_dir = calc_wm_ref_stats(
            wm_ref_speed=wm_ref_speed,
            wm_ref_direction=wm_ref_direction,
            epsg=meta.product.grid.epsg,
            bounds=meta.product.geometry.bbox_native,
        )
        meta.product.wind_normalization = WindNormalizationMetadata(
            backscatter_measurement='sigma0',
            backscatter_convention='intensity ratio',
            reference_direction=wm_ref_mean_dir,
            reference_model=URL['windNormReferenceModel'],
            reference_speed=wm_ref_mean_speed,
            reference_type='sigma0-ref',
        )
    else:
        meta.product.wind_normalization = None


def calc_geolocation_accuracy(
        swath_identifier: str,
        ei_tif: str,
        etad: bool,
        decimals: int = 2
) -> float | None:
    """
    Calculates the radial root-mean-square error, which is a target requirement of the CARD4L NRB specification
    (Item 4.3). For more information see: https://s1ard.readthedocs.io/en/latest/general/geoaccuracy.html.
    Currently, only the Copernicus DEM is supported.

    Parameters
    ----------
    swath_identifier
        Swath identifier dependent on acquisition mode.
    ei_tif
        Path to the annotation GeoTIFF layer 'Ellipsoidal Incident Angle' of the current product.
    etad
        Was the ETAD correction applied?
    decimals
        Number of decimal places to round the calculated rRMSE value to. Default is 2.

    Returns
    -------
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


def calc_pslr_islr(
        annotation_dict: dict[str, dict[str, etree.ElementTree]],
        decimals: int = 2
) -> tuple[float | None, float | None]:
    """
    Extracts all values for Peak Side Lobe Ratio (PSLR) and Integrated
    Side Lobe Ratio (ISLR) from the annotation metadata of a scene and
    calculates the mean value for all swaths.

    Parameters
    ----------
    annotation_dict
        A dictionary of Sentinel-1 annotation file XML objects
        (key `annotation` of the dictionary returned by :func:`get_src_meta`).
    decimals
        Number of decimal places to round the calculated values to. Default is 2.

    Returns
    -------
        a tuple with the following values. Each can be set to None if not found in
        the metadata.

        - pslr: Mean PSLR value for all swaths of the scene.
        - islr: Mean ISLR value for all swaths of the scene.
    """
    out = []
    for key in ['crossCorrelationPslr', 'crossCorrelationIslr']:
        slr_dict = find_in_annotation(annotation_dict=annotation_dict,
                                      pattern=f'.//{key}',
                                      out_type='float',
                                      per_pol=True)
        # Mean values per swath
        slr_mean = [np.nanmean(v) for inner in slr_dict.values()
                    for v in inner.values() if v is not None]
        
        # Mean value for all swaths
        slr = round(np.nanmean(slr_mean).item(), decimals) if len(slr_mean) > 0 else None
        out.append(slr)
    return tuple(out)


def calc_wm_ref_stats(
        wm_ref_speed: list[str],
        wm_ref_direction: list[str],
        epsg: int,
        bounds: list[int | float],
        resolution: int = 915
) -> tuple[float, float]:
    """
    Calculates the mean wind model reference speed and direction for the wind model annotation layer.

    Parameters
    ----------
    wm_ref_speed
        List of paths pointing to the wind model reference speed files.
    wm_ref_speed
        List of paths pointing to the wind model reference direction files.
    epsg
        The EPSG code of the current MGRS tile.
    bounds
        The bounds of the current MGRS tile.
    resolution
        The resolution of the wind model reference files in meters. Default is 915.

    Returns
    -------
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


def copy_src_meta(
        ard_dir: str,
        src_ids: list[ID]
) -> None:
    """
    Copies the original metadata of the source scenes to the ARD product
    directory.

    Parameters
    ----------
    ard_dir
        A path pointing to the current ARD product directory.
    src_ids
        List of :class:`~pyroSAR.drivers.ID` objects of all source scenes that overlap with the current MGRS tile.
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


def find_in_annotation(
        annotation_dict: dict[str, dict[str, etree.ElementTree]],
        pattern: str,
        single: bool = False,
        per_pol: bool = False,
        out_type: Literal["int", "float", "str"] = "str"
) -> (dict[str, dict[str, list[int | float | str] | int | float | str | None]] |
      list[int | float | str] | int | float | str | None):
    """
    Search for a pattern in all XML annotation files provided and return the results.

    Parameters
    ----------
    annotation_dict
        A dictionary of Sentinel-1 annotation file XML objects
        (key `annotation` of the dictionary returned by :func:`get_src_meta`).
    pattern
        The pattern to search for in each annotation file.
    single
        If True, the results found in each annotation file are expected to be
        the same and therefore only a single value will be returned instead of
        a dictionary. If the results differ, an error is raised.
    per_pol
        group the results per polarization? In this case, the dictionary becomes
        nested with the outer dict having the swaths as keys and the inner the
        polarizations.
    out_type
        Output type to convert the results to.

    Returns
    -------
        Either a dictionary with results per annotation file (swath, polarization)
        or a single result, which is identical across annotation files. None is
        set if the value is not found.
    """
    out = defaultdict(defaultdict)
    for swath in annotation_dict.keys():
        for pol in annotation_dict[swath].keys():
            sub = annotation_dict[swath][pol]
            swaths = [x.text for x in sub.findall('.//swathProcParams/swath')]
            items = sub.findall(pattern)
            if len(items) > 0:
                parent = items[0].getparent().tag
                if parent in ['azimuthProcessing', 'rangeProcessing']:
                    for i, val in enumerate(items):
                        out[swaths[i]][pol] = val.text
                else:
                    out[swath][pol] = [x.text for x in items]
                    if len(out[swath][pol]) == 1:
                        out[swath][pol] = out[swath][pol][0]
            else:
                out[swath][pol] = None
    out = defaultdict_to_dict(out)
    
    def _convert(
            obj: dict[str, list[str] | str] | list[str] | str,
            out_type: str
    ):
        if isinstance(obj, dict):
            return {k: _convert(v, out_type) for k, v in obj.items()}
        if isinstance(obj, list):
            return [_convert(x, out_type) for x in obj]
        elif isinstance(obj, str):
            if out_type == 'float':
                return float(obj)
            if out_type == 'int':
                return int(obj)
            else:
                raise ValueError(f"unknown output type: {out_type}")
        elif obj is None:
            return obj
        else:
            raise TypeError(f'got an unexpected type: {type(obj)}')
    
    if out_type != 'str':
        out = _convert(out, out_type)
    
    err_msg = 'Search result for pattern "{}" expected to be the same in all annotation files.'
    if single:
        vals = [x for inner in out.values() for x in inner.values()]
        all_same = all(x == vals[0] for x in vals)
        if not all_same:
            raise RuntimeError(err_msg.format(pattern))
        return vals[0]
    
    err_msg = 'Search result for pattern "{}" expected to be the same for all polarizations.'
    if not per_pol:
        vals = [[x for x in inner.values()] for inner in out.values()]
        all_same = all(all(x == sub[0] for x in sub) for sub in vals)
        if not all_same:
            raise RuntimeError(err_msg.format(pattern))
        out = {k: inner[next(iter(inner))] for k, inner in out.items()}
    
    return out


def get_osv_info(
        sid: ID
) -> tuple[str, str] | tuple[None, None]:
    """
    Get information about the used OSV file.
    First, this function attempts to find an auxiliary OSV file matching the scene.
    If found, its name is returned. If not, it is assumed that it is not yet available
    and processing was performed using the OSVs found in the source product.
    In this case, the metadata is searched for the name of an auxiliary OSV file
    used during L1 generation. If found, its name is returned.

    Parameters
    ----------
    sid
        The pyroSAR scene ID object

    Returns
    -------
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


def get_prod_meta(
        tif: str,
        src_ids: list[ID],
        sar_dir: str,
        processor_name: str
) -> dict[str, Any]:
    """
    Get ARD product metadata.
    Returns a metadata dictionary, which is generated from a measurement
    GeoTIFF file of the ARD product and function `get_metadata` of the SAR
    processor.

    Parameters
    ----------
    tif
        The path to a measurement GeoTIFF file of the ARD product.
    src_ids
        List of objects of all source products.
    sar_dir
        A path pointing to the processed SAR datasets of the product.
    processor_name
        The name of the SAR processor. Needed for reading processing metadata,
        e.g. :func:`s1ard.snap.get_metadata`.

    Returns
    -------
        A dictionary containing metadata for the product scene.
    """
    processor = load_processor(processor_name)
    out = dict()
    coord_list = [
        coordinate
        for sid in src_ids
        for coordinate in sid.meta["coordinates"]
    ]
    
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
        
        # Calculate the number of nodata border pixels based on the footprint
        # of the source product(s)
        with get_geometry(coordinates=coord_list, crs=4326) as srcvec:
            ras_srcvec = rasterize(vectorobject=srcvec, reference=ras,
                                   burn_values=[1])
            arr_srcvec = ras_srcvec.array()
            out['nodata_borderpx'] = int(np.count_nonzero(np.isnan(arr_srcvec)))
    
    proc_meta = processor.get_metadata(scene=src_ids[0].scene, outdir=sar_dir)
    out['ML_nRgLooks'] = proc_meta['rlks']
    out['ML_nAzLooks'] = proc_meta['azlks']
    return out


# ElementTree is a cython class, which may behave like a function and may
# thus not define __or__, which is needed for the | typing operator.
# Using Union instead.
def get_src_meta(
        sid: ID
) -> dict[Literal["manifest", "annotation"], Union[etree.ElementTree, dict[str, dict[str, etree.ElementTree]]]]:
    """
    Retrieve the manifest and annotation XML data of a scene as a dictionary
    using an :class:`pyroSAR.drivers.ID` object.

    Parameters
    ----------
    sid
        A pyroSAR :class:`~pyroSAR.drivers.ID` object generated with e.g.
        :func:`pyroSAR.drivers.identify`.

    Returns
    -------
        A dictionary containing the parsed `etree.ElementTree` objects for
        the manifest and annotation XML files.
    """
    annotation_files = sid.findfiles(r'^s1[abcd].*-[vh]{2}-.*\.xml$')
    
    annotation = defaultdict(defaultdict)
    for file in annotation_files:
        base = os.path.basename(file)
        pol = re.search('[vh]{2}', base).group()
        swath = re.search('-(iw[1-3]*|ew[1-5]*|s[1-6])', base).group(1)
        content = etree.fromstring(sid.getFileObj(file).getvalue())
        annotation[swath.upper()][pol.upper()] = content
    
    with sid.getFileObj(sid.findfiles('manifest.safe')[0]) as input_man:
        manifest = etree.fromstring(input_man.getvalue())
    
    return {'manifest': manifest,
            'annotation': defaultdict_to_dict(annotation)}


def meta_dict(
        config: dict[str, Any],
        prod_meta: dict[str, Any],
        src_ids: list[ID],
        compression: str
) -> ARDMetadata:
    """
    Creates an :class:`~cesard.metadata.model.ARDMetadata` object containing
    metadata for a product scene as well as its source scenes. The object can
    then be used by :func:`cesard.metadata.xml.parse` and
    :func:`cesard.metadata.stac.parse` to generate OGC XML and STAC JSON
    metadata files, respectively.

    Parameters
    ----------
    config
        Dictionary of the parsed config parameters for the current process.
    prod_meta
        a metadata dictionary as returned by :func:`s1ard.ard.product_info`
    src_ids
        List of :class:`~pyroSAR.drivers.ID` objects of all source scenes
        that overlap with the current MGRS tile.
    compression
        The compression type applied to raster files of the product.

    Returns
    -------
        An :class:`~cesard.metadata.model.ARDMetadata` object containing
        metadata for the product as well as its source scenes.
    """
    
    def _as_utc(value: datetime | str) -> datetime:
        if isinstance(value, str):
            value = dateparse(value)
        if value.tzinfo is None or value.utcoffset() is None:
            return value.replace(tzinfo=timezone.utc)
        return value.astimezone(timezone.utc)
    
    src_sid = {}
    src_xml = {}
    for sid in src_ids:
        uid = os.path.basename(sid.scene).split('.')[0][-4:]
        src_sid[uid] = sid
        src_xml[uid] = get_src_meta(sid=sid)
    sid0 = next(iter(src_sid.values()))
    
    ref_tif = finder(prod_meta['dir_ard_product'], ['[hv]{2}-[gs]-lin.tif$'], regex=True)[0]
    np_tifs = finder(prod_meta['dir_ard_product'], ['-np-[hv]{2}.tif$'], regex=True)
    ei_tif = finder(prod_meta['dir_ard_product'], ['-ei.tif$'], regex=True)
    prod_meta.update(get_prod_meta(
        tif=ref_tif,
        src_ids=src_ids,
        sar_dir=config['processing']['sar_dir'],
        processor_name=config['processing']['processor']
    ))
    op_mode = prod_meta['mode']
    
    # COMMON metadata
    pid_lookup = {'S1A': 'A', 'S1B': 'B', 'S1C': 'C', 'S1D': 'D'}
    platform_identifier = pid_lookup[sid0.sensor]
    platform_short_name = 'Sentinel-1'
    platform_full_name = f'{platform_short_name.lower()}{platform_identifier.lower()}'
    common = CommonMetadata(
        antenna_look_direction='RIGHT',
        constellation='sentinel-1',
        instrument_short_name='C-SAR',
        operational_mode=op_mode,
        orbit_direction={'A': 'ascending', 'D': 'descending'}[sid0.orbit],
        orbit_mean_altitude=693000.0,
        orbit_number_absolute=sid0.meta['orbitNumber_abs'],
        orbit_number_relative=sid0.meta['orbitNumber_rel'],
        platform_full_name=platform_full_name,
        platform_identifier=platform_identifier,
        platform_reference=URL['platformReference'][platform_full_name],
        platform_short_name=platform_short_name,
        polarizations=tuple(sid0.polarizations),
        polarization_mode=prod_meta['polarization'][0],
        processing_level='L1C',
        radar_band='C',
        radar_center_frequency=5405000000.0,
        sensor_type='RADAR',
        swath_identifier=op_mode,
        wrs_longitude_grid=sid0.meta['orbitNumbers_rel']['start'],
    )
    
    # PRODUCT metadata
    if prod_meta['product_type'] == 'ORB':
        card4l_link = URL['card4l_orb']
        card4l_version = '1.0'
    else:
        card4l_link = URL['card4l_nrb']
        card4l_version = '5.5'
    
    dem = DEM_MAP[config['processing']['dem_type']]
    dem_gsd = dem['gsd']
    if isinstance(dem_gsd, str):
        value, unit = re.search(
            '(^[0-9]+) (m|arcsec)$',
            dem_gsd
        ).groups()
        dem_gsd = GroundSamplingDistance(
            value=float(value),
            unit=unit,
        )
    else:
        raise ValueError('DEM GSD must be a string')
    
    if (len(ei_tif) == 1 and
            sid0.product == 'SLC' and
            'copernicus' in config['processing']['dem_type'].lower()):
        geo_corr_accuracy = calc_geolocation_accuracy(
            swath_identifier=op_mode,
            ei_tif=ei_tif[0],
            etad=False
        )
    else:
        geo_corr_accuracy = None
    
    if geo_corr_accuracy is not None:
        geo_corr_reference = URL['geoCorrAccuracyReference']
        geo_corr_rmse = geo_corr_accuracy
    else:
        geo_corr_reference = NOT_IMPLEMENTED_TEXT
        geo_corr_rmse = NOT_IMPLEMENTED_NUMBER
    
    processing_center = config['metadata']['processing_center']
    if processing_center is None:
        processing_center = NOT_IMPLEMENTED_TEXT
    
    rad_acc_reference = URL['radiometricAccuracyReference']
    if rad_acc_reference is None:
        rad_acc_reference = NOT_IMPLEMENTED_TEXT
    
    product_name = (
        'Ocean Radar Backscatter'
        if prod_meta['product_type'] == 'ORB'
        else 'Normalised Radar Backscatter'
    )
    noise_removal_applied = True
    
    product = ProductMetadata(
        name=product_name,
        product_type=prod_meta['product_type'],
        acquisition_type='NOMINAL',
        status='PLANNED',
        access=config['metadata']['access_url'],
        doi=config['metadata']['doi'],
        license=config['metadata']['licence'],
        time_start=_as_utc(prod_meta['start']),
        time_stop=_as_utc(prod_meta['stop']),
        time_created=_as_utc(prod_meta['proc_time']),
        geometry=GeometryMetadata(
            bbox=tuple(float(x) for x in prod_meta['geom']['bbox']),
            bbox_native=tuple(float(x) for x in prod_meta['geom']['bbox_native']),
            geometry=dict(prod_meta['geom']['geometry']),
        ),
        grid=GridMetadata(
            epsg=prod_meta['epsg'],
            wkt=prod_meta['wkt'],
            rows=prod_meta['rows'],
            columns=prod_meta['cols'],
            pixel_spacing_row=float(prod_meta['res'][1]),
            pixel_spacing_column=float(prod_meta['res'][0]),
            transform=tuple(float(x) for x in prod_meta['transform']),
            mgrs_id=prod_meta['tile'],
            pixel_coordinate_convention='upper-left',
            number_of_border_pixels=prod_meta['nodata_borderpx'],
            definition_reference=URL['grid_definition'],
            convention_reference=URL['grid_convention'],
        ),
        processing=ProcessingMetadata(
            facility=processing_center,
            date=_as_utc(prod_meta['proc_time']),
            mode='PROTOTYPE',
            processor='s1ard',
            software={
                str(key): str(value)
                for key, value in version_dict(
                    processor_name=config['processing']['processor']
                ).items()
            },
        ),
        card4l=Card4lMetadata(
            specification=prod_meta['product_type'],
            version=card4l_version,
            document=card4l_link,
        ),
        compression=CompressionMetadata(
            type=compression,
            z_errors=LERC_ERR_THRES,
        ),
        dem=DEMMetadata(
            name=config['processing']['dem_type'].replace(' II', ''),
            type=dem['type'],
            reference=dem['ref'],
            access=dem['access'],
            gsd=dem_gsd,
            resampling_method='bilinear',
            egm_reference=dem['egm'],
            egm_resampling_method='bilinear' if dem['egm'] is not None else None,
        ),
        backscatter=BackscatterMetadata(
            measurement='gamma0' if re.search('g-lin', ref_tif) else 'sigma0',
            convention='linear power',
            conversion_equation='10*log10(DN)',
            range_number_of_looks=round(prod_meta['ML_nRgLooks'], 2),
            azimuth_number_of_looks=round(prod_meta['ML_nAzLooks'], 2),
            equivalent_number_of_looks=calc_enl(tif=ref_tif),
        ),
        geometric_correction=GeometricCorrectionMetadata(
            algorithm=URL['geoCorrAlgorithm'],
            resampling_method='bilinear',
            accuracy=GeometricAccuracyMetadata(
                type='gtc',
                eastern=AxisAccuracy(
                    bias=NOT_IMPLEMENTED_NUMBER,
                    standard_deviation=NOT_IMPLEMENTED_NUMBER,
                ),
                northern=AxisAccuracy(
                    bias=NOT_IMPLEMENTED_NUMBER,
                    standard_deviation=NOT_IMPLEMENTED_NUMBER,
                ),
                radial_rmse=geo_corr_rmse,
                reference=geo_corr_reference,
            ),
        ),
        radiometric_accuracy=RadiometricAccuracyMetadata(
            absolute=NOT_IMPLEMENTED_NUMBER,
            relative=NOT_IMPLEMENTED_NUMBER,
            reference=rad_acc_reference,
        ),
        noise_removal=NoiseRemovalMetadata(
            applied=noise_removal_applied,
            algorithm=URL['noiseRemovalAlgorithm'] if noise_removal_applied else None,
        ),
        rtc_algorithm=URL['RTCAlgorithm'],
        number_of_acquisitions=len(src_sid),
        speckle_filter_applied=None,
        ellipsoidal_height=NOT_IMPLEMENTED_NUMBER,
        wind_normalization=None,
    )
    
    # SOURCE metadata
    sources = {}
    for uid, sid in src_sid.items():
        nsmap = src_xml[uid]['manifest'].nsmap
        
        swath_ids = find_in_annotation(
            annotation_dict=src_xml[uid]['annotation'],
            pattern='.//swathProcParams/swath'
        )
        swaths = []
        for item in swath_ids.values():
            if isinstance(item, list):
                swaths.extend(item)
            else:
                swaths.append(item)
        
        with sid.geometry() as vec:
            geom = geometry_from_vec(vectorobject=vec)
        
        res_mode = re.match(
            re.compile(sid.pattern), os.path.basename(sid.file)
        ).groupdict()['resolution']
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
        else:
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
            result = find_in_annotation(
                annotation_dict=src_xml[uid]['annotation'],
                pattern=pattern,
                out_type=out_type
            )
            results.append(result)
        (az_look_bandwidth, rg_look_bandwidth, az_num_looks, rg_num_looks,
         az_px_spacing, rg_px_spacing, inc) = results
        
        inc_vals = dissolve(list(inc.values()))
        lut_applied = find_in_annotation(
            annotation_dict=src_xml[uid]['annotation'],
            pattern='.//applicationLutId',
            single=True
        )
        pslr, islr = calc_pslr_islr(annotation_dict=src_xml[uid]['annotation'])
        
        def _read_manifest(pattern, attrib=None):
            obj = src_xml[uid]['manifest'].find(pattern, nsmap)
            return obj.attrib[attrib] if attrib is not None else obj.text
        
        asc_node_time = _as_utc(_read_manifest('.//s1:ascendingNodeTime'))
        proc_date = _as_utc(
            _read_manifest('.//safe:processing', attrib='stop')
        )
        proc_name = _read_manifest('.//safe:software', attrib='name')
        proc_version = _read_manifest('.//safe:software', attrib='version')
        fac_org = _read_manifest('.//safe:facility', attrib='organisation')
        fac_name = _read_manifest('.//safe:facility', attrib='name')
        
        if len(np_tifs) > 0:
            performance_estimates = calc_performance_estimates(files=np_tifs)
            noise_equivalent_intensity_type = 'sigma0'
        else:
            performance_estimates = {
                pol: {'minimum': None, 'mean': None, 'maximum': None}
                for pol in common.polarizations
            }
            noise_equivalent_intensity_type = None
        
        inc_min = round(min(inc_vals), 2)
        inc_max = round(max(inc_vals), 2)
        inc_mid = round(inc_max - ((inc_max - inc_min) / 2), 2)
        osv_base, osv_descr = get_osv_info(sid)
        
        sources[uid] = SourceMetadata(
            filename=sid.file,
            product_type=product_type,
            data_geometry=data_geometry,
            acquisition_type='NOMINAL',
            status='ARCHIVED',
            access=URL['source_access'],
            doi=URL['source_doi'],
            time_start=_as_utc(sid.start),
            time_stop=_as_utc(sid.stop),
            geometry=GeometryMetadata(
                bbox=tuple(float(x) for x in geom['bbox']),
                geometry=dict(geom['geometry']),
            ),
            processing=ProcessingMetadata(
                facility=f"{fac_org} {fac_name}".replace(' -', ''),
                date=proc_date,
                mode='NOMINAL',
                processor=proc_name,
                software={proc_name: proc_version},
            ),
            orbit=SourceOrbitMetadata(
                ascending_node_date=asc_node_time,
                start_time_from_ascending_node=float(
                    _read_manifest('.//s1:startTimeANX')
                ),
                completion_time_from_ascending_node=float(
                    _read_manifest('.//s1:stopTimeANX')
                ),
                major_cycle_id=sid.meta['cycleNumber'],
                datatake_id=int(_read_manifest('.//s1sarl1:missionDataTakeID')),
                data_access=URL['orbitDataAccess'],
                data_source=osv_descr,
                state_vector=osv_base,
            ),
            swaths=tuple(swaths),
            azimuth=SwathAxisMetadata(
                look_bandwidth=az_look_bandwidth,
                number_of_looks=az_num_looks,
                pixel_spacing=az_px_spacing,
                resolution=res_az,
            ),
            range=SwathAxisMetadata(
                look_bandwidth=rg_look_bandwidth,
                number_of_looks=rg_num_looks,
                pixel_spacing=rg_px_spacing,
                resolution=res_rg,
            ),
            incidence_angle=IncidenceAngleMetadata(
                minimum=inc_min,
                maximum=inc_max,
                mid_swath=inc_mid,
            ),
            instrument_azimuth_angle=round(sid.meta['heading'], 2),
            lut_applied=lut_applied,
            sensor_calibration=URL['sensorCalibration'],
            polarimetric_calibration_matrices=None,
            faraday_mean_rotation_angle=None,
            faraday_rotation_reference=URL['faradayRotationReference'],
            ionosphere_indicator=None,
            performance=SourcePerformanceMetadata(
                noise_equivalent_intensity_type=noise_equivalent_intensity_type,
                estimates={
                    pol: PerformanceEstimate(**values)
                    for pol, values in performance_estimates.items()
                },
                equivalent_number_of_looks=enl,
                integrated_side_lobe_ratio=islr,
                peak_side_lobe_ratio=pslr,
            ),
        )
    
    return ARDMetadata(
        common=common,
        product=product,
        sources=sources,
    )
