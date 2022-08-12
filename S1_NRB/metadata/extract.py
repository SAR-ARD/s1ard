import os
import re
import math
from lxml import etree
from datetime import datetime
import numpy as np
from statistics import median
from pyroSAR.snap.auxil import parse_recipe
from spatialist import Raster
from spatialist.ancillary import finder, dissolve
from spatialist.vector import wkt2vector, bbox
from spatialist.raster import rasterize
from osgeo import gdal
import S1_NRB
from S1_NRB.metadata.mapping import NRB_PATTERN, ITEM_MAP, RES_MAP, ORB_MAP, DEM_MAP, SLC_ACC_MAP

gdal.UseExceptions()


def get_prod_meta(product_id, tif, src_ids, snap_outdir):
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
    snap_outdir: str
        A path pointing to the SNAP processed datasets of the product.
    
    Returns
    -------
    dict
        A dictionary containing metadata for the product scene.
    """
    out = re.match(re.compile(NRB_PATTERN), product_id).groupdict()
    coord_list = [sid.meta['coordinates'] for sid in src_ids]
    
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
            out['extent_4326'] = vec.extent
            
            # Calculate number of nodata border pixels based on source scene(s) footprint
            ras_srcvec = rasterize(vectorobject=srcvec, reference=ras, burn_values=[1])
            arr_srcvec = ras_srcvec.array()
            out['nodata_borderpx'] = np.count_nonzero(np.isnan(arr_srcvec))
    
    pat = 'S1[AB]__(IW|EW|S[1-6])___(A|D)_[0-9]{8}T[0-9]{6}.+ML.+xml$'
    wf_path = finder(snap_outdir, [pat], regex=True)
    if len(wf_path) > 0:
        wf = parse_recipe(wf_path[0])
        rlks = int(wf['Multilook'].parameters['nRgLooks'])
        azlks = int(wf['Multilook'].parameters['nAzLooks'])
    else:
        rlks = azlks = 1
    
    src_xml = etree_from_sid(sid=src_ids[0])
    az_num_looks = find_in_annotation(annotation_dict=src_xml['annotation'],
                                      pattern='.//azimuthProcessing/numberOfLooks',
                                      out_type='int')
    rg_num_looks = find_in_annotation(annotation_dict=src_xml['annotation'],
                                      pattern='.//rangeProcessing/numberOfLooks',
                                      out_type='int')
    out['ML_nRgLooks'] = rlks * median(rg_num_looks.values())
    out['ML_nAzLooks'] = azlks * median(az_num_looks.values())
    return out


def vec_from_srccoords(coord_list):
    """
    Creates a single :class:`~spatialist.vector.Vector` object from a list of footprint coordinates of source scenes.
    
    Parameters
    ----------
    coord_list: list[list[tuple(float, float)]]
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


def etree_from_sid(sid):
    """
    Retrieve the manifest and annotation XML data of a scene as a dictionary using an :class:`~pyroSAR.drivers.ID`
    object.
    
    Parameters
    ----------
    sid:  :class:`pyroSAR.drivers.ID`
        A pyroSAR :class:`~pyroSAR.drivers.ID` object generated with :func:`pyroSAR.drivers.identify`.
    
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


def convert_coordinates(coords, stac=False):
    """
    Converts footprint coordinates that have been retrieved from the metadata of source SLC scenes stored in an
    :class:`~pyroSAR.drivers.ID` object OR a product extent retrieved using :class:`spatialist.vector.Vector.extent` to
    either `envelop` and `center` for usage in the XML metadata files or `bbox` and `geometry` for usage in STAC
    metadata files. The latter is returned if the optional parameter `stac` is set to True, else the former is returned.
    
    Parameters
    ----------
    coords: list[tuple(float, float)] or dict
        List of coordinate tuple pairs as retrieved from an :class:`~pyroSAR.drivers.ID` objects of source SLC scenes
        OR the product extent retrieved using :class:`spatialist.vector.Vector.extent` in the form of a dictionary with
        keys: xmin, xmax, ymin, ymax
    stac: bool, optional
        If set to True, `bbox` and `geometry` are returned for usage in STAC metadata file. If set to False (default)
        `envelop` and `center` are returned for usage in XML metadata files.
    
    Returns
    -------
    envelop: str
        Acquisition footprint coordinates for the XML element 'eop:Footprint/multiExtentOf'.
    center: str
        Acquisition center coordinates for the XML element 'eop:Footprint/centerOf'.
    
    Notes
    -------
    If `stac=True` the following results are returned instead of `envelop` and `center`:
    
    bbox: list[float]
        Acquisition bounding box for usage in STAC Items. Formatted in accordance with RFC 7946, section 5:
        https://datatracker.ietf.org/doc/html/rfc7946#section-5
    geometry: dict
        Acquisition footprint geometry for usage in STAC Items. Formatted in accordance with RFC 7946, section 3.1.:
        https://datatracker.ietf.org/doc/html/rfc7946#section-3.1
    """
    if isinstance(coords, (list, tuple)) and len(coords) == 4:
        c = coords
        x = [c[0][0], c[1][0], c[2][0], c[3][0]]
        y = [c[0][1], c[1][1], c[2][1], c[3][1]]
        xmin = min(x)
        xmax = max(x)
        ymin = min(y)
        ymax = max(y)
    elif isinstance(coords, dict) and len(coords.keys()) == 4:
        xmin = coords['xmin']
        xmax = coords['xmax']
        ymin = coords['ymin']
        ymax = coords['ymax']
        x = [xmin, xmin, xmax, xmax]
        y = [ymin, ymax, ymax, ymin]
    else:
        raise RuntimeError('Coordinates must be provided as a list of coordinate tuples OR as a dictionary with '
                           'keys xmin, xmax, ymin, ymax')
    if stac:
        bbox = [xmin, ymin, xmax, ymax]
        geometry = {'type': 'Polygon', 'coordinates': (((x[0], y[0]), (x[1], y[1]), (x[2], y[2]), (x[3], y[3]),
                                                        (x[0], y[0])),)}
        return bbox, geometry
    else:
        x_c = (xmax + xmin) / 2
        y_c = (ymax + ymin) / 2
        center = '{} {}'.format(y_c, x_c)
        envelop = '{} {} {} {} {} {} {} {} {} {}'.format(y[0], x[0], y[1], x[1], y[2], x[2], y[3], x[3], y[0], x[0])
        return center, envelop


def find_in_annotation(annotation_dict, pattern, single=False, out_type='str'):
    """
    Search for a pattern in all XML annotation files provided and return a dictionary of results.
    
    Parameters
    ----------
    annotation_dict: dict
        A dict of annotation files in the form: {'swath ID': `lxml.etree._Element` object}
    pattern: str
        The pattern to search for in each annotation file.
    single: bool, optional
        If True, the results found in each annotation file are expected to be the same and therefore only a single
        value will be returned instead of a dict. If the results differ, an error is raised. Default is False.
    out_type: str, optional
        Output type to convert the results to. Can be one of the following:
        
        - str (default)
        - float
        - int
    
    Returns
    -------
    out: dict
        A dictionary of the results containing a list for each of the annotation files. E.g.,
        {'swath ID': list[str, float or int]}
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
    
    def convert(obj, type):
        if isinstance(obj, list):
            return [convert(x, type) for x in obj]
        elif isinstance(obj, str):
            if type == 'float':
                return float(obj)
            if type == 'int':
                return int(obj)
    
    if out_type != 'str':
        for k, v in list(out.items()):
            out[k] = convert(v, out_type)
    
    err_msg = 'Search result for pattern "{}" expected to be the same in all annotation files.'
    if single:
        val = list(out.values())[0]
        for k in out:
            if out[k] != val:
                raise RuntimeError(err_msg.format(pattern))
        if out_type != 'str':
            return convert(val, out_type)
        else:
            return val
    else:
        return out


def calc_performance_estimates(files, ref_tif):
    """
    Calculates the performance estimates specified in CARD4L NRB 1.6.9 for all noise power images if available.
    
    Parameters
    ----------
    files: list[str]
        List of paths pointing to the noise power images the estimates should be calculated for.
    ref_tif: str
        A path pointing to a reference product raster, which is used to get spatial information about the current MGRS
        tile.
    
    Returns
    -------
    out: dict
        Dictionary containing the calculated estimates for each available polarization.
    """
    out = {}
    with Raster(ref_tif) as ref:
        for f in files:
            pol = re.search('[VH]{2}', f).group().upper()
            with bbox(ref.extent, crs=ref.epsg) as vec:
                with Raster(f)[vec] as ras:
                    arr = ras.array()
                    # The following need to be of type float, not numpy.float32 in order to be JSON serializable
                    _min = float(np.nanmin(arr))
                    _max = float(np.nanmax(arr))
                    _mean = float(np.nanmean(arr))
                    del arr
            out[pol] = {'minimum': _min,
                        'maximum': _max,
                        'mean': _mean}
    return out


def extract_pslr_islr(annotation_dict):
    """
    Extracts all values for Peak Side Lobe Ratio (PSLR) and Integrated Side Lobe Ratio (ISLR) from the annotation
    metadata of a scene and calculates the mean value for all swaths.
    
    Parameters
    ----------
    annotation_dict: dict
        A dictionary of annotation files in the form: {'swath ID':`lxml.etree._Element` object}
    
    Returns
    -------
    pslr: float
        Mean PSLR value for all swaths of the scene.
    islr: float
        Mean ISLR value for all swaths of the scene.
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
        A path to a GeoTIFF file of the currently processed NRB product.
    
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


def calc_geolocation_accuracy(swath_identifier, ei_tif, dem_type, etad):
    """
    Calculates the radial root mean square error, which is a target requirement of the CARD4L NRB specification
    (Item 4.3). For more information see: https://s1-nrb.readthedocs.io/en/latest/general/geoaccuracy.html
    
    Parameters
    ----------
    swath_identifier: str
        Swath identifier dependent on acquisition mode.
    ei_tif: str
        Path to the annotation GeoTIFF layer 'Ellipsoidal Incident Angle' of the current product.
    dem_type: str
        The DEM type used for processing.
    etad: bool
        Was the ETAD correction applied?
    
    Returns
    -------
    rmse_planar: float
        The calculated rRMSE value rounded to two decimal places.
    """
    if 'copernicus' not in dem_type.lower():
        return None
    
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


def meta_dict(config, target, src_ids, snap_datasets, proc_time, start, stop, compression):
    """
    Creates a dictionary containing metadata for a product scene, as well as its source scenes. The dictionary can then
    be utilized by :func:`~S1_NRB.metadata.xml.parse` and :func:`~S1_NRB.metadata.stac.parse` to generate XML and STAC
    JSON metadata files, respectively.
    
    Parameters
    ----------
    config: dict
        Dictionary of the parsed config parameters for the current process.
    target: str
        A path pointing to the NRB product scene being created.
    src_ids: list[pyroSAR.drivers.ID]
        List of :class:`~pyroSAR.drivers.ID` objects of all source scenes that overlap with the current MGRS tile.
    snap_datasets: list[str]
        List of output files processed with :func:`pyroSAR.snap.util.geocode` that match the source SLC scenes
        overlapping with the current MGRS tile.
    proc_time: datetime.datetime
        The processing time object used to generate the unique product identifier.
    start: datetime.datetime
        The product start time.
    stop: datetime.datetime
        The product stop time.
    compression: str
        The compression type applied to raster files of the product.
    
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
        src_xml[uid] = etree_from_sid(sid=sid)
    sid0 = src_sid[list(src_sid.keys())[0]]  # first key/first file; used to extract some common metadata
    swath_id = re.search('_(IW|EW|S[1-6])_', os.path.basename(sid0.file)).group().replace('_', '')
    
    ref_tif = finder(target, ['[hv]{2}-g-lin.tif$'], regex=True)[0]
    np_tifs = finder(target, ['-np-[hv]{2}.tif$'], regex=True)
    ei_tif = finder(target, ['-ei.tif$'], regex=True)[0]
    
    product_id = os.path.basename(target)
    prod_meta = get_prod_meta(product_id=product_id, tif=ref_tif, src_ids=src_ids,
                              snap_outdir=os.path.dirname(snap_datasets[0]))
    
    xml_center, xml_envelop = convert_coordinates(coords=prod_meta['extent_4326'])
    stac_bbox, stac_geometry = convert_coordinates(coords=prod_meta['extent_4326'], stac=True)
    stac_bbox_native = convert_coordinates(coords=prod_meta['extent'], stac=True)[0]
    
    dem_type = config['dem_type']
    dem_access = DEM_MAP[dem_type]['access']
    dem_ref = DEM_MAP[dem_type]['ref']
    dem_subtype = DEM_MAP[dem_type]['type']
    egm_ref = DEM_MAP[dem_type]['egm']
    dem_name = dem_type.replace(' II', '')
    
    tups = [(ITEM_MAP[key]['suffix'], ITEM_MAP[key]['z_error']) for key in ITEM_MAP.keys()]
    z_err_dict = dict(tups)
    
    geocorr_acc = calc_geolocation_accuracy(swath_identifier=swath_id, ei_tif=ei_tif, dem_type=dem_type,
                                            etad=config['etad'])
    
    # Common metadata (sorted alphabetically)
    meta['common']['antennaLookDirection'] = 'RIGHT'
    meta['common']['constellation'] = 'sentinel-1'
    meta['common']['instrumentShortName'] = 'C-SAR'
    meta['common']['operationalMode'] = prod_meta['mode']
    meta['common']['orbitDirection'] = {'A': 'ascending', 'D': 'descending'}[sid0.orbit]
    meta['common']['orbitMeanAltitude'] = '{:.2e}'.format(693000)
    meta['common']['orbitNumber'] = str(sid0.meta['orbitNumbers_abs']['stop'])
    meta['common']['orbitNumbers_abs'] = sid0.meta['orbitNumbers_abs']
    meta['common']['orbitNumbers_rel'] = sid0.meta['orbitNumbers_rel']
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
    
    # Product metadata (sorted alphabetically)
    meta['prod']['access'] = config['meta']['access_url']
    meta['prod']['ancillaryData_KML'] = 'https://sentinel.esa.int/documents/247904/1955685/S2A_OPER_GIP_TILPAR_MPC__' \
                                        '20151209T095117_V20150622T000000_21000101T000000_B00.kml'
    meta['prod']['acquisitionType'] = 'NOMINAL'
    meta['prod']['azimuthNumberOfLooks'] = prod_meta['ML_nAzLooks']
    meta['prod']['backscatterConvention'] = 'linear power'
    meta['prod']['backscatterConversionEq'] = '10*log10(DN)'
    meta['prod']['backscatterMeasurement'] = 'gamma0'
    meta['prod']['card4l-link'] = 'https://ceos.org/ard/files/PFS/NRB/v5.5/CARD4L-PFS_NRB_v5.5.pdf'
    meta['prod']['card4l-version'] = '5.5'
    meta['prod']['crsEPSG'] = str(prod_meta['epsg'])
    meta['prod']['crsWKT'] = prod_meta['wkt']
    meta['prod']['compression_type'] = compression
    meta['prod']['compression_zerrors'] = z_err_dict
    meta['prod']['demEGMReference'] = egm_ref
    meta['prod']['demEGMResamplingMethod'] = 'bilinear'
    meta['prod']['demName'] = dem_name
    meta['prod']['demReference'] = dem_ref
    meta['prod']['demResamplingMethod'] = 'bilinear'
    meta['prod']['demType'] = dem_subtype
    meta['prod']['demAccess'] = dem_access
    meta['prod']['doi'] = config['meta']['doi']
    meta['prod']['ellipsoidalHeight'] = None
    meta['prod']['fileBitsPerSample'] = '32'
    meta['prod']['fileByteOrder'] = 'little-endian'
    meta['prod']['fileDataType'] = 'float'
    meta['prod']['fileFormat'] = 'COG'
    meta['prod']['speckleFilterApplied'] = False
    meta['prod']['geoCorrAccuracyEasternBias'] = None
    meta['prod']['geoCorrAccuracyEasternSTDev'] = None
    meta['prod']['geoCorrAccuracyNorthernBias'] = None
    meta['prod']['geoCorrAccuracyNorthernSTDev'] = None
    meta['prod']['geoCorrAccuracy_rRMSE'] = geocorr_acc
    meta['prod']['geoCorrAccuracyReference'] = 'https://s1-nrb.readthedocs.io/en/v{}/general/geoaccuracy.html' \
                                               ''.format(S1_NRB.__version__)
    meta['prod']['geoCorrAccuracyType'] = 'slant-range'
    meta['prod']['geoCorrAlgorithm'] = 'https://sentinel.esa.int/documents/247904/1653442/' \
                                       'Guide-to-Sentinel-1-Geocoding.pdf'
    meta['prod']['geoCorrResamplingMethod'] = 'bilinear'
    meta['prod']['geom_stac_bbox_native'] = stac_bbox_native
    meta['prod']['geom_stac_bbox_4326'] = stac_bbox
    meta['prod']['geom_stac_geometry_4326'] = stac_geometry
    meta['prod']['geom_xml_center'] = xml_center
    meta['prod']['geom_xml_envelope'] = xml_envelop
    meta['prod']['griddingConventionURL'] = 'http://www.mgrs-data.org/data/documents/nga_mgrs_doc.pdf'
    meta['prod']['griddingConvention'] = 'Military Grid Reference System (MGRS)'
    meta['prod']['licence'] = config['meta']['licence']
    meta['prod']['mgrsID'] = prod_meta['mgrsID']
    meta['prod']['NRApplied'] = True if len(np_tifs) > 0 else False
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
    meta['prod']['productName'] = 'Normalised Radar Backscatter'
    meta['prod']['productName-short'] = 'NRB'
    meta['prod']['pxSpacingColumn'] = str(prod_meta['res'][0])
    meta['prod']['pxSpacingRow'] = str(prod_meta['res'][1])
    meta['prod']['radiometricAccuracyAbsolute'] = None
    meta['prod']['radiometricAccuracyRelative'] = None
    meta['prod']['radiometricAccuracyReference'] = None
    meta['prod']['rangeNumberOfLooks'] = prod_meta['ML_nRgLooks']
    meta['prod']['RTCAlgorithm'] = 'https://doi.org/10.1109/Tgrs.2011.2120616'
    meta['prod']['status'] = 'PLANNED'
    meta['prod']['timeCreated'] = proc_time
    meta['prod']['timeStart'] = start
    meta['prod']['timeStop'] = stop
    meta['prod']['transform'] = prod_meta['transform']
    
    # Source metadata
    for uid in list(src_sid.keys()):
        nsmap = src_xml[uid]['manifest'].nsmap
        
        swath_ids = find_in_annotation(annotation_dict=src_xml[uid]['annotation'],
                                       pattern='.//swathProcParams/swath')
        swaths = []
        for item in swath_ids.values():
            if isinstance(item, list):
                swaths.extend(item)
            else:
                swaths.append(item)
        osv = src_sid[uid].getOSV(returnMatch=True, osvType=['POE', 'RES'], useLocal=True)
        
        coords = src_sid[uid].meta['coordinates']
        xml_center, xml_envelop = convert_coordinates(coords=coords)
        stac_bbox, stac_geometry = convert_coordinates(coords=coords, stac=True)
        
        az_look_bandwidth = find_in_annotation(annotation_dict=src_xml[uid]['annotation'],
                                               pattern='.//azimuthProcessing/lookBandwidth',
                                               out_type='float')
        az_num_looks = find_in_annotation(annotation_dict=src_xml[uid]['annotation'],
                                          pattern='.//azimuthProcessing/numberOfLooks')
        az_px_spacing = find_in_annotation(annotation_dict=src_xml[uid]['annotation'],
                                           pattern='.//azimuthPixelSpacing',
                                           out_type='float')
        inc = find_in_annotation(annotation_dict=src_xml[uid]['annotation'],
                                 pattern='.//geolocationGridPoint/incidenceAngle',
                                 out_type='float')
        inc_vals = dissolve(inc.values())
        lut_applied = find_in_annotation(annotation_dict=src_xml[uid]['annotation'],
                                         pattern='.//applicationLutId', single=True)
        pslr, islr = extract_pslr_islr(annotation_dict=src_xml[uid]['annotation'])
        np_files = [ds for ds in snap_datasets if re.search('_NE[BGS]Z', ds) is not None]
        rg_look_bandwidth = find_in_annotation(annotation_dict=src_xml[uid]['annotation'],
                                               pattern='.//rangeProcessing/lookBandwidth',
                                               out_type='float')
        rg_num_looks = find_in_annotation(annotation_dict=src_xml[uid]['annotation'],
                                          pattern='.//rangeProcessing/numberOfLooks')
        rg_px_spacing = find_in_annotation(annotation_dict=src_xml[uid]['annotation'],
                                           pattern='.//rangePixelSpacing',
                                           out_type='float')
        
        def read_manifest(pattern, attrib=None):
            obj = src_xml[uid]['manifest'].find(pattern, nsmap)
            if attrib is not None:
                return obj.attrib[attrib]
            else:
                return obj.text
        
        # (sorted alphabetically)
        meta['source'][uid] = {}
        meta['source'][uid]['access'] = 'https://scihub.copernicus.eu'
        meta['source'][uid]['acquisitionType'] = 'NOMINAL'
        meta['source'][uid]['ascendingNodeDate'] = read_manifest('.//s1:ascendingNodeTime')
        meta['source'][uid]['azimuthLookBandwidth'] = az_look_bandwidth
        meta['source'][uid]['azimuthNumberOfLooks'] = az_num_looks
        meta['source'][uid]['azimuthPixelSpacing'] = az_px_spacing
        op_mode = meta['common']['operationalMode']
        if re.search('S[1-6]', op_mode):
            res_az = {op_mode: RES_MAP['SM']['azimuthResolution'][op_mode]}
            res_rg = {op_mode: RES_MAP['SM']['rangeResolution'][op_mode]}
        else:
            res_az = RES_MAP[op_mode]['azimuthResolution']
            res_rg = RES_MAP[op_mode]['rangeResolution']
        meta['source'][uid]['azimuthResolution'] = res_az
        if src_sid[uid].meta['product'] == 'GRD':
            meta['source'][uid]['dataGeometry'] = 'ground range'
        else:
            meta['source'][uid]['dataGeometry'] = 'slant range'
        meta['source'][uid]['datatakeID'] = read_manifest('.//s1sarl1:missionDataTakeID')
        meta['source'][uid]['doi'] = 'https://sentinel.esa.int/documents/247904/1877131/' \
                                     'Sentinel-1-Product-Specification'
        meta['source'][uid]['faradayMeanRotationAngle'] = None
        meta['source'][uid]['faradayRotationReference'] = None
        meta['source'][uid]['filename'] = src_sid[uid].file
        meta['source'][uid]['geom_stac_bbox_4326'] = stac_bbox
        meta['source'][uid]['geom_stac_geometry_4326'] = stac_geometry
        meta['source'][uid]['geom_xml_center'] = xml_center
        meta['source'][uid]['geom_xml_envelop'] = xml_envelop
        meta['source'][uid]['incidenceAngleMax'] = np.max(inc_vals)
        meta['source'][uid]['incidenceAngleMin'] = np.min(inc_vals)
        meta['source'][uid]['incidenceAngleMidSwath'] = np.max(inc_vals) - ((np.max(inc_vals) - np.min(inc_vals)) / 2)
        meta['source'][uid]['instrumentAzimuthAngle'] = str(src_sid[uid].meta['heading'])
        meta['source'][uid]['ionosphereIndicator'] = None
        meta['source'][uid]['lutApplied'] = lut_applied
        meta['source'][uid]['majorCycleID'] = str(src_sid[uid].meta['cycleNumber'])
        meta['source'][uid]['orbitStateVector'] = os.path.basename(osv).replace('.zip', '')
        for orb in list(ORB_MAP.keys()):
            if orb in meta['source'][uid]['orbitStateVector']:
                meta['source'][uid]['orbitDataSource'] = ORB_MAP[orb]
        meta['source'][uid]['orbitDataAccess'] = 'https://scihub.copernicus.eu/gnss'
        if len(np_files) > 0:
            meta['source'][uid]['perfEstimates'] = calc_performance_estimates(files=np_files, ref_tif=ref_tif)
            meta['source'][uid]['perfNoiseEquivalentIntensityType'] = 'sigma0'
        else:
            stats = {stat: None for stat in ['minimum', 'mean', 'maximum']}
            pe = {pol: stats for pol in meta['common']['polarisationChannels']}
            meta['source'][uid]['perfEstimates'] = pe
            meta['source'][uid]['perfNoiseEquivalentIntensityType'] = None
        meta['source'][uid]['perfEquivalentNumberOfLooks'] = 1
        meta['source'][uid]['perfIntegratedSideLobeRatio'] = islr
        meta['source'][uid]['perfPeakSideLobeRatio'] = pslr
        meta['source'][uid]['polCalMatrices'] = None
        fac_org = read_manifest('.//safe:facility', attrib='organisation')
        fac_name = read_manifest('.//safe:facility', attrib='name')
        meta['source'][uid]['processingCenter'] = f"{fac_org} {fac_name}".replace(' -', '')
        meta['source'][uid]['processingDate'] = read_manifest('.//safe:processing', attrib='stop')
        meta['source'][uid]['processingLevel'] = read_manifest('.//safe:processing', attrib='name')
        meta['source'][uid]['processorName'] = read_manifest('.//safe:software', attrib='name')
        meta['source'][uid]['processorVersion'] = read_manifest('.//safe:software', attrib='version')
        meta['source'][uid]['processingMode'] = 'NOMINAL'
        meta['source'][uid]['productType'] = src_sid[uid].meta['product']
        meta['source'][uid]['rangeLookBandwidth'] = rg_look_bandwidth
        meta['source'][uid]['rangeNumberOfLooks'] = rg_num_looks
        meta['source'][uid]['rangePixelSpacing'] = rg_px_spacing
        meta['source'][uid]['azimuthResolution'] = res_az
        meta['source'][uid]['rangeResolution'] = res_rg
        meta['source'][uid]['sensorCalibration'] = 'https://sentinel.esa.int/web/sentinel/technical-guides/' \
                                                   'sentinel-1-sar/sar-instrument/calibration'
        meta['source'][uid]['status'] = 'ARCHIVED'
        meta['source'][uid]['swaths'] = swaths
        meta['source'][uid]['timeCompletionFromAscendingNode'] = str(float(read_manifest('.//s1:stopTimeANX')))
        meta['source'][uid]['timeStartFromAscendingNode'] = str(float(read_manifest('.//s1:startTimeANX')))
        meta['source'][uid]['timeStart'] = datetime.strptime(src_sid[uid].start, '%Y%m%dT%H%M%S')
        meta['source'][uid]['timeStop'] = datetime.strptime(src_sid[uid].stop, '%Y%m%dT%H%M%S')
    
    return meta
