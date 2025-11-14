import os
import re
from pathlib import Path
from dateutil.parser import parse as dateparse
from packaging.version import Version
from pystac_client import Client
from pystac_client.stac_api_io import StacApiIO
from spatialist.vector import Vector
from shapely.geometry import shape
from pyroSAR import identify_many
from cesard.ancillary import date_to_utc, buffer_time
from cesard.search import asf_select
import logging

log = logging.getLogger('s1ard')


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
        lookup_platform = {
            'S1A': 'sentinel-1a',
            'S1B': 'sentinel-1b',
            'S1C': 'sentinel-1c',
            'S1D': 'sentinel-1Dd'
        }
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
