import os
import re
import inspect
from pathlib import Path
from dateutil.parser import parse as dateparse
from packaging.version import Version
from datetime import timedelta
from pystac_client import Client
from pystac_client.stac_api_io import StacApiIO
from spatialist.vector import Vector, crsConvert, wkt2vector
from shapely.geometry import shape
import asf_search as asf
from pyroSAR import identify_many, ID
from s1ard.ancillary import date_to_utc, buffer_time, combine_polygons
from s1ard.tile_extraction import aoi_from_tile, tile_from_aoi
import logging

log = logging.getLogger('s1ard')


class ASF(ID):
    """
    Simple SAR metadata handler for scenes in the ASF archive. The interface is consistent with the driver classes in
    :mod:`pyroSAR.drivers` but does not implement the full functionality due to limited content of the CMR
    metadata catalog. Registered attributes:
    
    - acquisition_mode
    - coordinates
    - frameNumber
    - orbit
    - orbitNumber_abs
    - orbitNumber_rel
    - polarizations
    - product
    - projection
    - sensor
    - start
    - stop
    """
    
    def __init__(self, meta):
        self.scene = meta['properties']['url']
        self._meta = meta
        self.meta = self.scanMetadata()
        super(ASF, self).__init__(self.meta)
    
    def __lt__(self, other):
        if not isinstance(other, ASF):
            return NotImplemented
        return self.outname_base() < other.outname_base()
    
    def scanMetadata(self):
        meta = dict()
        meta['acquisition_mode'] = self._meta['properties']['beamModeType']
        meta['coordinates'] = [tuple(x) for x in self._meta['geometry']['coordinates'][0]]
        fname = os.path.splitext(self._meta['properties']['fileName'])[0]
        meta['frameNumber'] = fname[-11:-5]
        meta['orbit'] = self._meta['properties']['flightDirection'][0]
        meta['orbitNumber_abs'] = self._meta['properties']['orbit']
        meta['orbitNumber_rel'] = self._meta['properties']['pathNumber']
        meta['polarizations'] = self._meta['properties']['polarization'].split('+')
        product = self._meta['properties']['processingLevel']
        meta['product'] = re.search('GRD|SLC|SM|OCN', product).group()
        meta['projection'] = crsConvert(4326, 'wkt')
        meta['sensor'] = self._meta['properties']['platform'].replace('entinel-', '')
        start = self._meta['properties']['startTime']
        stop = self._meta['properties']['stopTime']
        pattern = '%Y%m%dT%H%M%S'
        meta['start'] = dateparse(start).strftime(pattern)
        meta['stop'] = dateparse(stop).strftime(pattern)
        meta['spacing'] = None
        meta['samples'] = None
        meta['lines'] = None
        meta['cycleNumber'] = None
        meta['sliceNumber'] = None
        meta['totalSlices'] = None
        return meta


class STACArchive(object):
    """
    Search for scenes in a SpatioTemporal Asset Catalog.
    Scenes are expected to be unpacked with a folder suffix .SAFE.
    The interface is kept consistent with :class:`~s1ard.search.ASFArchive`,
    :class:`~s1ard.search.STACParquetArchive` and :class:`pyroSAR.drivers.Archive`.
    
    Parameters
    ----------
    url: str
        the catalog URL
    collections: str or list[str]
        the catalog collection(s) to be searched
    timeout: int
        the allowed timeout in seconds
    max_retries: int or None
        the number of times to retry requests. Set to None to disable retries.
    
    See Also
    --------
    pystac_client.Client.open
    pystac_client.stac_api_io.StacApiIO
    """
    
    def __init__(self, url, collections, timeout=60, max_retries=20):
        self.url = url
        self.timeout = timeout
        self.max_tries = max_retries
        self._open_catalog()
        if isinstance(collections, str):
            self.collections = [collections]
        elif isinstance(collections, list):
            self.collections = collections
        else:
            raise TypeError("'collections' must be of type str or list")
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
    
    def _filter_duplicates(self, values):
        tmp = sorted(values, key=lambda x: os.path.basename(x[-2]))
        pattern = '([0-9A-Z_]{16})_([0-9T]{15})_([0-9T]{15})'
        keep = []
        i = 0
        while i < len(tmp):
            group = [tmp[i]]
            base = os.path.basename(tmp[i][-2])
            match1 = re.search(pattern, base).groups()
            j = i + 1
            while j < len(tmp):
                base = os.path.basename(tmp[j][-2])
                match2 = re.search(pattern, base).groups()
                if match1 == match2:
                    group.append(tmp[j])
                    j += 1
                else:
                    break
            if len(group) > 1:
                tproc = [x[-1] for x in group]
                keep.append(group[tproc.index(max(tproc))])
            else:
                keep.append(group[0])
            i = j
        return keep
    
    def _open_catalog(self):
        stac_api_io = StacApiIO(max_retries=self.max_tries)
        self.catalog = Client.open(url=self.url,
                                   stac_io=stac_api_io,
                                   timeout=self.timeout)
    
    def close(self):
        del self.catalog
    
    def select(self, sensor=None, product=None, acquisition_mode=None,
               mindate=None, maxdate=None, frameNumber=None,
               vectorobject=None, date_strict=True, check_exist=True,
               return_value="scene"):
        """
        Select scenes from the catalog. Duplicates (same acquisition time) are filtered
        by returning only the last processed product. Used STAC property keys:

        - platform
        - start_datetime
        - end_datetime
        - created
        - sar:instrument_mode
        - sar:product_type
        - s1:datatake (custom)

        Parameters
        ----------
        sensor: str or list[str] or None
            S1A or S1B
        product: str or list[str] or None
            GRD or SLC
        acquisition_mode: str or list[str] or None
            IW, EW or SM
        mindate: str or datetime.datetime or None
            the minimum acquisition date; timezone-unaware dates are interpreted as UTC.
        maxdate: str or datetime.datetime or None
            the maximum acquisition date; timezone-unaware dates are interpreted as UTC.
        frameNumber: int or str or list[int or str] or None
            the data take ID in decimal (int) or hexadecimal (str) representation.
            Requires custom STAC key `s1:datatake`.
        vectorobject: spatialist.vector.Vector or None
            a geometry with which the scenes need to overlap. The object may only contain one feature.
        date_strict: bool
            treat dates as strict limits or also allow flexible limits to incorporate scenes
            whose acquisition period overlaps with the defined limit?

            - strict: start >= mindate & stop <= maxdate
            - not strict: stop >= mindate & start <= maxdate
        check_exist: bool
            check whether found files exist locally?
        return_value: str or List[str]
            the query return value(s). Options:
            
            - acquisition_mode: the sensor's acquisition mode, e.g., IW, EW, SM
            - frameNumber: the frame or datatake number
            - geometry_wkb: the scene's footprint geometry formatted as WKB
            - geometry_wkt: the scene's footprint geometry formatted as WKT
            - mindate: the acquisition start datetime in UTC formatted as YYYYmmddTHHMMSS
            - maxdate: the acquisition end datetime in UTC formatted as YYYYmmddTHHMMSS
            - product: the product type, e.g., SLC, GRD
            - scene: the scene's storage location path (default)
            - sensor: the satellite platform, e.g., S1A or S1B

        Returns
        -------
        list or list[tuple]
            If a single return_value is specified: list of values
            If multiple return_values are specified: list of tuples containing the requested values

        See Also
        --------
        pystac_client.Client.search
        """
        pars = locals()
        del pars['date_strict']
        del pars['check_exist']
        del pars['return_value']
        del pars['self']
        
        if isinstance(return_value, str):
            return_values = [return_value]
        else:
            return_values = return_value
        
        lookup = {'product': 'sar:product_type',
                  'acquisition_mode': 'sar:instrument_mode',
                  'mindate': 'start_datetime',
                  'maxdate': 'end_datetime',
                  'sensor': 'platform',
                  'frameNumber': 's1:datatake'}
        lookup_platform = {'S1A': 'sentinel-1a',
                           'S1B': 'sentinel-1b'}
        lookup_platform_reverse = {value: key for key, value
                                   in lookup_platform.items()}
        
        args = {'datetime': [None, None]}
        flt = {'op': 'and', 'args': []}
        for key in pars.keys():
            val = pars[key]
            if val is None:
                continue
            if key in ['mindate', 'maxdate']:
                val = date_to_utc(val, str_format='%Y-%m-%dT%H:%M:%SZ')
            if key == 'mindate':
                args['datetime'][0] = val
                if date_strict:
                    arg = {'op': '>=', 'args': [{'property': 'start_datetime'}, val]}
                else:
                    arg = {'op': '>=', 'args': [{'property': 'end_datetime'}, val]}
                flt['args'].append(arg)
            elif key == 'maxdate':
                args['datetime'][1] = val
                if date_strict:
                    arg = {'op': '<=', 'args': [{'property': 'end_datetime'}, val]}
                else:
                    arg = {'op': '<=', 'args': [{'property': 'start_datetime'}, val]}
                flt['args'].append(arg)
            elif key == 'vectorobject':
                if isinstance(val, Vector):
                    if val.nfeatures > 1:
                        raise RuntimeError("'vectorobject' may only contain one feature")
                    with val.clone() as vec:
                        vec.reproject(4326)
                        feat = vec.getFeatureByIndex(0)
                        json = feat.ExportToJson(as_object=True)['geometry']
                        feat = None
                        args['intersects'] = json
                else:
                    raise TypeError('argument vectorobject must be of type spatialist.vector.Vector')
            else:
                if isinstance(val, (str, int)):
                    val = [val]
                val_format = []
                for v in val:
                    if key == 'sensor':
                        v = lookup_platform[v]
                    if key == 'frameNumber' and isinstance(v, int):
                        v = '{:06X}'.format(v)  # convert to hexadecimal
                    val_format.append(v)
                if len(val_format) == 1:
                    arg = {'op': '=', 'args': [{'property': lookup[key]}, val_format[0]]}
                else:
                    arg = {'op': 'in', 'args': [{'property': lookup[key]}, val_format]}
                flt['args'].append(arg)
        if len(flt['args']) == 0:
            flt = None
        if args['datetime'] == [None, None]:
            args['datetime'] = None
        result = self.catalog.search(collections=self.collections,
                                     filter=flt, max_items=None,
                                     **args)
        result = list(result.items())
        out = []
        for item in result:
            assets = item.assets
            ref = assets[list(assets.keys())[0]]
            href = ref.href
            scene = href[:re.search(r'\.SAFE', href).end()]
            scene = re.sub('^file://', '', scene)
            if Path(scene).exists():
                scene = os.path.realpath(scene)
            else:
                if check_exist:
                    raise RuntimeError('scene does not exist locally:', scene)
            
            # prepare return values
            values = []
            for key in return_values:
                if key == "scene":
                    values.append(scene)
                elif key in lookup.keys():
                    value = item.properties[lookup[key]]
                    if key in ['mindate', 'maxdate']:
                        value = dateparse(value).strftime('%Y%m%dT%H%M%S')
                    if key == 'sensor':
                        value = lookup_platform_reverse[value]
                    values.append(value)
                # reverse coordinate order to be consistent with ASF return
                elif key == "geometry_wkt":
                    values.append(shape(item.geometry).wkt)
                elif key == "geometry_wkb":
                    values.append(shape(item.geometry).wkb)
                else:
                    raise ValueError(f"Invalid return value: {key}")
            
            values.append(scene)
            created = dateparse(item.properties['created'])
            values.append(created)
            
            out.append(values)
        
        out = self._filter_duplicates(out)
        
        def reduce(i):
            o = i[:-2]
            return o[0] if len(o) == 1 else tuple(o)
        
        out = [reduce(x) for x in out]
        return out


class STACParquetArchive(object):
    """
    Search for scenes in STAC geoparquet dump.
    Scenes are expected to be unpacked with a folder suffix .SAFE.
    The interface is kept consistent with :class:`~s1ard.search.ASFArchive`,
    :class:`~s1ard.search.STACArchive` and :class:`pyroSAR.drivers.Archive`.

    Parameters
    ----------
    files: str
        the file search pattern, e.g. `/path/to/*parquet`
    """
    
    def __init__(self, files):
        self.files = files
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        return
    
    def close(self):
        pass
    
    def select(self, sensor=None, product=None, acquisition_mode=None,
               mindate=None, maxdate=None, frameNumber=None,
               vectorobject=None, date_strict=True, return_value='scene'):
        """
        Select scenes from a STAC catalog's geoparquet dump. Used STAC property keys:

        - platform
        - start_datetime
        - end_datetime
        - sar:instrument_mode
        - sar:product_type
        - s1:datatake (custom)

        Parameters
        ----------
        sensor: str or list[str] or None
            S1A or S1B
        product: str or list[str] or None
            GRD or SLC
        acquisition_mode: str or list[str] or None
            IW, EW or SM
        mindate: str or datetime.datetime or None
            the minimum acquisition date; timezone-unaware dates are interpreted as UTC.
        maxdate: str or datetime.datetime or None
            the maximum acquisition date; timezone-unaware dates are interpreted as UTC.
        frameNumber: int or str or list[int or str] or None
            the data take ID in decimal (int) or hexadecimal (str) representation.
            Requires custom STAC key `s1:datatake`.
        vectorobject: spatialist.vector.Vector or None
            a geometry with which the scenes need to overlap
        date_strict: bool
            treat dates as strict limits or also allow flexible limits to incorporate scenes
            whose acquisition period overlaps with the defined limit?

            - strict: start >= mindate & stop <= maxdate
            - not strict: stop >= mindate & start <= maxdate
        return_value: str or List[str]
            the query return value(s). Options:
            
            - acquisition_mode: the sensor's acquisition mode, e.g., IW, EW, SM
            - frameNumber: the frame or datatake number
            - geometry_wkb: the scene's footprint geometry formatted as WKB
            - geometry_wkt: the scene's footprint geometry formatted as WKT
            - mindate: the acquisition start datetime in UTC formatted as YYYYmmddTHHMMSS
            - maxdate: the acquisition end datetime in UTC formatted as YYYYmmddTHHMMSS
            - product: the product type, e.g., SLC, GRD
            - scene: the scene's storage location path (default)
            - sensor: the satellite platform, e.g., S1A or S1B

        Returns
        -------
        List[str] or List[tuple[str]]
            the selected return value(s). Depending on whether a single or multiple
            values have been defined for `return_value`, the returned list will
            contain strings or tuples.
        
        See Also
        --------
        stac_geoparquet.arrow.to_parquet
        duckdb.query
        """
        pars = locals()
        try:
            import duckdb
        except ImportError:
            raise ImportError("this method requires 'duckdb>=1.1.1' to be installed")
        ddb_version = Version(duckdb.__version__)
        ddb_version_req = Version('1.1.1')
        if ddb_version < ddb_version_req:
            raise ImportError("duckdb version must be >= 1.1.1")
        
        duckdb.install_extension('spatial')
        duckdb.load_extension('spatial')
        
        del pars['self']
        del pars['date_strict']
        del pars['return_value']
        
        lookup = {'product': 'sar:product_type',
                  'acquisition_mode': 'sar:instrument_mode',
                  'mindate': 'start_datetime',
                  'maxdate': 'end_datetime',
                  'sensor': 'platform',
                  'frameNumber': 's1:datatake'}
        lookup_platform = {'S1A': 'sentinel-1a',
                           'S1B': 'sentinel-1b',
                           'S1C': 'sentinel-1c',
                           'S1D': 'sentinel-1d'}
        return_value_mapping = {
            "geometry_wkb": "ST_AsWKB(geometry)",
            "geometry_wkt": "ST_AsText(geometry)",
            "mindate": f"STRFTIME({lookup['mindate']} AT TIME ZONE 'UTC', '%Y%m%dT%H%M%S')",
            "maxdate": f"STRFTIME({lookup['maxdate']} AT TIME ZONE 'UTC', '%Y%m%dT%H%M%S')",
            "scene": "replace(json_extract_string(assets::json, '$.folder.href'), 'file://', '')"
        }
        for k, v in lookup.items():
            if k not in return_value_mapping:
                return_value_mapping[k] = f"\"{lookup[k]}\""
        
        return_values = return_value if isinstance(return_value, list) else [return_value]
        for return_value in return_values:
            if return_value not in return_value_mapping:
                raise ValueError(f"unsupported return value '{return_value}'.\n"
                                 f"supported options: {list(return_value_mapping.keys())}")
        
        terms = []
        for key in pars.keys():
            val = pars[key]
            if val is None:
                continue
            if key in ['mindate', 'maxdate']:
                val = date_to_utc(val, str_format='%Y-%m-%dT%H:%M:%SZ')
            if key == 'mindate':
                if date_strict:
                    terms.append(f'"start_datetime" >= \'{val}\'')
                else:
                    terms.append(f'"end_datetime" >= \'{val}\'')
            elif key == 'maxdate':
                if date_strict:
                    terms.append(f'"end_datetime" <= \'{val}\'')
                else:
                    terms.append(f'"start_datetime" <= \'{val}\'')
            elif key == 'vectorobject':
                if isinstance(val, Vector):
                    if val.nfeatures > 1:
                        raise RuntimeError("'vectorobject' may only contain one feature")
                    with val.clone() as tmp:
                        tmp.reproject(4326)
                        wkt = tmp.convert2wkt(set3D=False)[0]
                    terms.append(f'ST_Intersects(geometry, '
                                 f'ST_GeomFromText(\'{wkt}\'))')
                else:
                    raise TypeError('argument vectorobject must be of type spatialist.vector.Vector')
            else:
                if isinstance(val, (str, int)):
                    val = [val]
                val_format = []
                for v in val:
                    if key == 'sensor':
                        v = lookup_platform[v]
                    if key == 'frameNumber' and isinstance(v, int):
                        v = '{:06X}'.format(v)  # convert to hexadecimal
                    val_format.append(v)
                if len(val_format) == 1:
                    subterm = f'"{lookup[key]}" = \'{val_format[0]}\''
                else:
                    subterm = f'"{lookup[key]}" IN {tuple(val_format)}'
                terms.append(subterm)
        sql_where = ' AND '.join(terms)
        sql_return_value = ', '.join([return_value_mapping[x] for x in return_values])
        sql_query = f"""
        SELECT {sql_return_value}
        FROM '{self.files}' WHERE {sql_where}
        """
        result = duckdb.query(sql_query).fetchall()
        if 'sensor' in return_values:
            lookup_platform_reverse = {value: key for key, value in lookup_platform.items()}
            sensor_index = return_values.index('sensor')
            for i, item in enumerate(result):
                item_new = list(item)
                item_new[sensor_index] = lookup_platform_reverse[item[sensor_index]]
                result[i] = tuple(item_new)
        if len(return_values) == 1:
            out = [x[0] for x in result]
        else:
            out = result
        return sorted(out)


class ASFArchive(object):
    """
    Search for scenes in the Alaska Satellite Facility (ASF) catalog.
    The interface is kept consistent with :class:`~s1ard.search.STACArchive`,
    :class:`~s1ard.search.STACParquetArchive` and :class:`pyroSAR.drivers.Archive`.
    """
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        return
    
    @staticmethod
    def select(sensor=None, product=None, acquisition_mode=None, mindate=None,
               maxdate=None, vectorobject=None, date_strict=True, return_value='scene'):
        """
        Select scenes from the ASF catalog. This is a simple wrapper around the function
        :func:`~s1ard.search.asf_select` to be consistent with the interfaces of the
        other search classes.

        Parameters
        ----------
        sensor: str or list[str] or None
            S1A or S1B
        product: str or list[str] or None
            GRD or SLC
        acquisition_mode: str or list[str] or None
            IW, EW or SM
        mindate: str or datetime.datetime or None
            the minimum acquisition date; timezone-unaware dates are interpreted as UTC.
        maxdate: str or datetime.datetime or None
            the maximum acquisition date; timezone-unaware dates are interpreted as UTC.
        vectorobject: spatialist.vector.Vector or None
            a geometry with which the scenes need to overlap. The object may only contain one feature.
        date_strict: bool
            treat dates as strict limits or also allow flexible limits to incorporate scenes
            whose acquisition period overlaps with the defined limit?
            
            - strict: start >= mindate & stop <= maxdate
            - not strict: stop >= mindate & start <= maxdate
        return_value: str or list[str]
            the metadata return value; see :func:`~s1ard.search.asf_select` for details
        
        See Also
        --------
        asf_select
        
        Returns
        -------
        list[str or tuple[str] or ASF]
            the scene metadata attributes as specified with `return_value`;
            see :func:`~s1ard.search.asf_select` for details
        """
        return asf_select(sensor, product, acquisition_mode, mindate, maxdate, vectorobject,
                          return_value=return_value, date_strict=date_strict)


def asf_select(sensor=None, product=None, acquisition_mode=None, mindate=None,
               maxdate=None, vectorobject=None, date_strict=True, return_value='scene'):
    """
    Search scenes in the Alaska Satellite Facility (ASF) data catalog.
    This is a simple interface to the
    `asf_search <https://github.com/asfadmin/Discovery-asf_search>`_ package.
    
    Parameters
    ----------
    sensor: str or None
        S1A|S1B|S1C|S1D
    product: str or None
        GRD|SLC
    acquisition_mode: str or None
        IW|EW|SM
    mindate: str or datetime.datetime or None
        the minimum acquisition date; timezone-unaware dates are interpreted as UTC.
    maxdate: str or datetime.datetime or None
        the maximum acquisition date; timezone-unaware dates are interpreted as UTC.
    vectorobject: spatialist.vector.Vector or None
        a geometry with which the scenes need to overlap. The object may only contain one feature.
    date_strict: bool
        treat dates as strict limits or also allow flexible limits to incorporate scenes
        whose acquisition period overlaps with the defined limit?
        
        - strict: start >= mindate & stop <= maxdate
        - not strict: stop >= mindate & start <= maxdate
    
    return_value: str or List[str]
        the query return value(s). Options:
        
        - acquisition_mode: the sensor's acquisition mode, e.g., IW, EW, SM
        - frameNumber: the frame or datatake number
        - geometry_wkb: the scene's footprint geometry formatted as WKB
        - geometry_wkt: the scene's footprint geometry formatted as WKT
        - mindate: the acquisition start datetime in UTC formatted as YYYYmmddTHHMMSS
        - maxdate: the acquisition end datetime in UTC formatted as YYYYmmddTHHMMSS
        - product: the product type, e.g., SLC, GRD
        - scene: the scene's storage location path (default)
        - sensor: the satellite platform, e.g., S1A or S1B

    Returns
    -------
    list[str or tuple[str] or ASF]
        the scene metadata attributes as specified with `return_value`; the return type
        is a list of strings, tuples, or :class:`~s1ard.search.ASF` objects depending on
        whether `return_type` is of type string, list or :class:`~s1ard.search.ASF`.
    
    """
    if isinstance(return_value, str):
        return_values = [return_value]
    else:
        return_values = return_value
    
    if product == 'GRD':
        processing_level = ['GRD_HD', 'GRD_MD', 'GRD_MS', 'GRD_HS', 'GRD_FD']
    else:
        processing_level = product
    if acquisition_mode == 'SM':
        beam_mode = ['S1', 'S2', 'S3', 'S4', 'S5', 'S6']
    else:
        beam_mode = acquisition_mode
    if vectorobject is not None:
        if vectorobject.nfeatures > 1:
            raise RuntimeError("'vectorobject' contains more than one feature.")
        with vectorobject.clone() as geom:
            geom.reproject(4326)
            geometry = geom.convert2wkt(set3D=False)[0]
    else:
        geometry = None
    
    start = date_to_utc(mindate, as_datetime=True)
    stop = date_to_utc(maxdate, as_datetime=True)
    
    lookup_platform = {'S1A': 'Sentinel-1A',
                       'S1B': 'Sentinel-1B',
                       'S1C': 'Sentinel-1C',
                       'S1D': 'Sentinel-1D'}
    platform = lookup_platform[sensor] if sensor is not None else None
    
    result = asf.search(platform=platform,
                        processingLevel=processing_level,
                        beamMode=beam_mode,
                        start=start,
                        end=stop,
                        intersectsWith=geometry).geojson()
    features = result['features']
    
    def date_extract(item, key):
        return date_to_utc(date=item['properties'][key], as_datetime=True)
    
    if date_strict:
        features = [x for x in features
                    if start <= date_extract(x, 'startTime')
                    and date_extract(x, 'stopTime') <= stop]
    
    features = sorted([ASF(x) for x in features])
    
    out = []
    for item in features:
        values = []
        for key in return_values:
            if key == 'ASF':
                values.append(item)
            elif key == 'mindate':
                values.append(getattr(item, 'start'))
            elif key == 'maxdate':
                values.append(getattr(item, 'stop'))
            elif key == 'geometry_wkb':
                with item.geometry() as vec:
                    value = vec.to_geopandas().to_wkb()['geometry'][0]
                    values.append(value)
            elif key == 'geometry_wkt':
                with item.geometry() as vec:
                    value = vec.to_geopandas().to_wkt()['geometry'][0]
                    values.append(value)
            elif hasattr(item, key):
                values.append(getattr(item, key))
            else:
                raise ValueError(f'invalid return value: {key}')
        if len(return_values) == 1:
            out.append(values[0])
        else:
            out.append(tuple(values))
    return out


def scene_select(archive, aoi_tiles=None, aoi_geometry=None, return_value='scene', **kwargs):
    """
    Central scene search utility. Selects scenes from a database and returns their file names
    together with the MGRS tile names for which to process ARD products.
    The list of MGRS tile names is either identical to the list provided with `aoi_tiles`,
    the list of all tiles overlapping with `aoi_geometry` or `vectorobject` (via `kwargs`),
    or the list of all tiles overlapping with an initial scene search result if no geometry
    has been defined via `aoi_tiles` or `aoi_geometry`. In the latter (most complex) case,
    the search procedure is as follows:
    
     - perform a first search matching all other search parameters
     - derive all MGRS tile geometries overlapping with the selection
     - derive the minimum and maximum acquisition times of the selection as search parameters
       `mindate` and `maxdate`
     - extend the `mindate` and `maxdate` search parameters by one minute
     - perform a second search with the extended time range and the derived MGRS tile geometries
     - filter the search result to scenes overlapping with the initial time range (if defined
       via `mindate` or `maxdate`)
    
    As consequence, if one defines the search parameters to only return one scene, the neighboring
    acquisitions will also be returned. This is because the scene overlaps with a set of MGRS
    tiles of which many or all will also overlap with these neighboring acquisitions. To ensure
    full coverage of all MGRS tiles, the neighbors of the scene in focus have to be processed too.
    
    This function has three ways to define search geometries. In order of priority overriding others:
    `aoi_tiles` > `aoi_geometry` > `vectorobject` (via `kwargs`). In the latter two cases, the search
    geometry is extended to the common footprint of all MGRS tiles overlapping with the initial geometry
    to ensure full coverage of all tiles.
    
    Parameters
    ----------
    archive: pyroSAR.drivers.Archive or STACArchive or STACParquetArchive or ASFArchive
        an open scene archive connection
    aoi_tiles: list[str] or None
        a list of MGRS tile names for spatial search
    aoi_geometry: str or None
        the name of a vector geometry file for spatial search
    return_value: str or List[str]
        the query return value(s). Default 'scene': return the scene's storage location path.
        See the documentation of `archive.select` for options.
    kwargs
        further search arguments passed to the `select` method of `archive`.
        The `date_strict` argument has no effect. Whether an ARD product is strictly in the defined
        time range cannot be determined by this function, and it thus has to add a time buffer.
        When `date_strict=True`, more scenes will be filtered out in the last step described above.

    Returns
    -------
    tuple[list[str or tuple[str]], list[str]]
    
     - the list of return values; single value:string, multiple values: tuple
     - the list of MGRS tiles
    
    """
    args = kwargs.copy()
    if 'mindate' in args.keys():
        args['mindate'] = date_to_utc(args['mindate'], as_datetime=True)
        mindate_init = args['mindate']
    else:
        mindate_init = None
    if 'maxdate' in args.keys():
        args['maxdate'] = date_to_utc(args['maxdate'], as_datetime=True)
        maxdate_init = args['maxdate']
    else:
        maxdate_init = None
    for key in ['acquisition_mode']:
        if key not in args.keys():
            args[key] = None
    
    if args['acquisition_mode'] == 'SM':
        args['acquisition_mode'] = ('S1', 'S2', 'S3', 'S4', 'S5', 'S6')
    
    signature = inspect.signature(archive.select)
    if 'return_value' not in signature.parameters:
        raise RuntimeError("the 'select' method of 'archive' does not take "
                           "a 'return_value' parameter")
    
    return_values = return_value if isinstance(return_value, list) else [return_value]
    if isinstance(archive, ASFArchive):
        args['return_value'] = 'ASF'
    else:
        args['return_value'] = return_value
    
    vec = None
    if aoi_tiles is not None:
        log.debug("reading geometries of 'aoi_tiles'")
        vec = aoi_from_tile(tile=aoi_tiles)
    elif aoi_geometry is not None:
        log.debug("extracting tiles overlapping with 'aoi_geometry'")
        with Vector(aoi_geometry) as geom:
            vec = tile_from_aoi(vector=geom,
                                return_geometries=True)
    elif 'vectorobject' in args.keys() and args['vectorobject'] is not None:
        log.debug("extracting tiles overlapping with 'vectorobject'")
        vec = tile_from_aoi(vector=args['vectorobject'],
                            return_geometries=True)
    if vec is not None:
        if not isinstance(vec, list):
            vec = [vec]
        
        if aoi_tiles is None:
            aoi_tiles = [x.mgrs for x in vec]
        log.debug(f"got {len(aoi_tiles)} tiles")
    
    # derive geometries and tiles from scene footprints
    if vec is None:
        log.debug("performing initial scene search without geometry constraint")
        args['return_value'] = ['mindate', 'maxdate', 'geometry_wkt']
        selection_tmp = archive.select(**args)
        log.debug(f'got {len(selection_tmp)} scenes')
        mindates, maxdates, geometries_init = zip(*selection_tmp)
        # The geometry of scenes crossing the antimeridian is stored as multipolygon.
        # Since the processor is currently not able to process these scenes, they are
        # removed in this step.
        geometries = [x for x in geometries_init if x.startswith('POLYGON')]
        if len(geometries) < len(geometries_init):
            log.debug(f'removed {len(geometries_init) - len(geometries)} '
                      f'scenes crossing the antimeridian')
        del selection_tmp
        
        log.debug(f"loading geometries")
        scenes_geom = [wkt2vector(x, srs=4326) for x in geometries]
        # select all tiles overlapping with the scenes for further processing
        log.debug("extracting all tiles overlapping with initial scene selection")
        vec = tile_from_aoi(vector=scenes_geom,
                            return_geometries=True)
        if not isinstance(vec, list):
            vec = [vec]
        aoi_tiles = [x.mgrs for x in vec]
        log.debug(f"got {len(aoi_tiles)} tiles")
        del scenes_geom
        
        args['mindate'] = min([date_to_utc(x, as_datetime=True) for x in mindates])
        args['maxdate'] = max([date_to_utc(x, as_datetime=True) for x in maxdates])
        del mindates, maxdates, geometries
    
    # extend the time range to fully cover all tiles
    # (one additional scene needed before and after each data take group)
    if 'mindate' in args.keys():
        args['mindate'] -= timedelta(minutes=1)
    if 'maxdate' in args.keys():
        args['maxdate'] += timedelta(minutes=1)
    
    args['return_value'] = return_values.copy()
    for key in ['mindate', 'maxdate']:
        if key not in args['return_value']:
            args['return_value'].append(key)
    
    log.debug("performing main scene search")
    with combine_polygons(vec, multipolygon=True) as combined:
        args['vectorobject'] = combined
        selection = archive.select(**args)
    del vec
    
    # reduce the selection to the time range defined by the user
    if mindate_init is not None or maxdate_init is not None:
        i = 0
        while i < len(selection):
            values = dict(zip(args['return_value'], selection[i]))
            start = date_to_utc(values['mindate'], as_datetime=True)
            stop = date_to_utc(values['maxdate'], as_datetime=True)
            delete = False
            if mindate_init is not None and stop < mindate_init:
                delete = True
            if maxdate_init is not None and start > maxdate_init:
                delete = True
            if delete:
                del selection[i]
            else:
                i += 1
    
    # sort the return values by the scene's basename
    rv_scene_key = args['return_value'].index('scene')
    selection = sorted(selection, key=lambda k: os.path.basename(k[rv_scene_key]))
    
    # reduce the return values to those defined by the user
    indices = [i for i, key in enumerate(args['return_value']) if key in return_values]
    if len(indices) == 1:
        selection = [scene[indices[0]] for scene in selection]
    else:
        selection = [tuple(scene[i] for i in indices) for scene in selection]
    log.debug(f"got {len(selection)} scenes")
    return selection, aoi_tiles


def collect_neighbors(archive, scene, stac_check_exist=True):
    """
    Collect a scene's neighboring acquisitions in a data take.
    
    Parameters
    ----------
    archive: pyroSAR.drivers.Archive or STACArchive or STACParquetArchive or ASFArchive
        an open scene archive connection
    scene: pyroSAR.drivers.ID
        the Sentinel-1 scene to be checked
    stac_check_exist: bool
        if `archive` is of type :class:`STACArchive`, check the local existence of the scenes?

    Returns
    -------
    list[str]
        the filenames/URLs of the neighboring scenes
    """
    start, stop = buffer_time(scene.start, scene.stop, seconds=2)
    
    kwargs = {'mindate': start, 'maxdate': stop, 'date_strict': False,
              'sensor': scene.sensor, 'product': scene.product,
              'acquisition_mode': scene.acquisition_mode,
              'return_value': ['scene', 'mindate', 'maxdate']}
    if isinstance(archive, STACArchive):
        kwargs['check_exist'] = stac_check_exist
    
    selection = archive.select(**kwargs)
    neighbors = [x for x in selection if x[1] != scene.start]
    if len(neighbors) > 2:
        # more than two neighbors can exist if multiple versions of the
        # datatake with different slicing exist.
        start_ref = dateparse(scene.start)
        stop_ref = dateparse(scene.stop)
        start_diff = []
        stop_diff = []
        for neighbor, start, stop in neighbors:
            start_diff.append(abs(start_ref - dateparse(stop)))
            stop_diff.append(abs(stop_ref - dateparse(start)))
        predecessor = neighbors[start_diff.index(min(start_diff))]
        successor = neighbors[stop_diff.index(min(stop_diff))]
        neighbors = [predecessor, successor]
    return [x[0] for x in neighbors]


def check_acquisition_completeness(archive, scenes):
    """
    Check presence of neighboring acquisitions.
    Check that for each scene a predecessor and successor can be queried
    from the database unless the scene is at the start or end of the data take.
    This ensures that no scene that could be covering an area of interest is missed
    during processing. In case a scene is suspected to be missing, the Alaska Satellite Facility (ASF)
    online catalog is cross-checked.
    An error will only be raised if the locally missing scene is present in the ASF catalog.
    It may happen that a neighbor is missing and this error is not raised if the scene is also
    missing on ASF.

    Parameters
    ----------
    archive: pyroSAR.drivers.Archive or STACArchive
        an open scene archive connection
    scenes: list[pyroSAR.drivers.ID]
        a list of scenes

    Returns
    -------

    Raises
    ------
    RuntimeError

    See Also
    --------
    s1ard.search.asf_select
    """
    messages = []
    for scene in scenes:
        log.debug(f'checking acquisition completeness for scene {scene.scene}')
        slice = scene.meta['sliceNumber']
        n_slices = scene.meta['totalSlices']
        groupsize = 3
        has_successor = True
        has_predecessor = True
        
        start, stop = buffer_time(scene.start, scene.stop, seconds=2)
        ref = None
        if slice == 0 or n_slices == 0:
            # NRT slicing mode
            ref = asf_select(sensor=scene.sensor,
                             product=scene.product,
                             acquisition_mode=scene.acquisition_mode,
                             mindate=start,
                             maxdate=stop,
                             return_value='scene')
            if len(ref) > 0:
                ref = [os.path.basename(x).replace('.zip', '.SAFE') for x in ref]
                match = [re.search(scene.pattern, x).groupdict() for x in ref]
                ref_start_min = min([x['start'] for x in match])
                ref_stop_max = max([x['stop'] for x in match])
                if ref_start_min == scene.start:
                    groupsize -= 1
                    has_predecessor = False
                if ref_stop_max == scene.stop:
                    groupsize -= 1
                    has_successor = False
            else:
                # don't assume neighbors if no products could be found of ASF
                has_successor = has_predecessor = False
        else:
            if slice == 1:  # first slice in the data take
                groupsize -= 1
                has_predecessor = False
            if slice == n_slices:  # last slice in the data take
                groupsize -= 1
                has_successor = False
        # Do another database selection to get the scene in question as well as its potential
        # predecessor and successor by adding an acquisition time buffer of two seconds.
        group = archive.select(sensor=scene.sensor,
                               product=scene.product,
                               acquisition_mode=scene.acquisition_mode,
                               mindate=start,
                               maxdate=stop,
                               date_strict=False)
        group = identify_many(group)
        # if the number of selected scenes is lower than the expected group size,
        # check whether the predecessor, the successor or both are missing by
        # cross-checking with the ASF database.
        if len(group) < groupsize:
            if ref is None:
                ref = asf_select(sensor=scene.sensor,
                                 product=scene.product,
                                 acquisition_mode=scene.acquisition_mode,
                                 mindate=start,
                                 maxdate=stop,
                                 return_value='scene')
            if len(ref) > 0:
                ref = [os.path.basename(x).replace('.zip', '.SAFE') for x in ref]
                match = [re.search(scene.pattern, x).groupdict() for x in ref]
                ref_start_min = min([x['start'] for x in match])
                ref_stop_max = max([x['stop'] for x in match])
                start_min = min([x.start for x in group])
                stop_max = max([x.stop for x in group])
                missing = []
                if ref_start_min < start < start_min and has_predecessor:
                    missing.append('predecessor')
                if stop_max < stop < ref_stop_max and has_successor:
                    missing.append('successor')
                if len(missing) > 0:
                    base = os.path.basename(scene.scene)
                    messages.append(f'{" and ".join(missing)} acquisition for scene {base}')
    if len(messages) != 0:
        text = '\n - '.join(messages)
        raise RuntimeError(f'missing the following scenes:\n - {text}')
