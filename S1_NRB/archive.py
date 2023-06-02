import os
import re
from lxml import etree
from pathlib import Path
from datetime import datetime
from pystac_client import Client
import pystac_client.exceptions
from spatialist import Vector


class STACArchive(object):
    """
    Search for scenes in a SpatioTemporal Asset Catalog.
    Scenes are expected to be unpacked with a folder suffix .SAFE.
    The interface is kept consistent with :class:`pyroSAR.drivers.Archive`.

    Parameters
    ----------
    url: str
        the catalog URL
    collections: str or list[str]
        the catalog collection(s) to be searched
    """
    
    def __init__(self, url, collections):
        self.max_tries = 3
        t = 1
        while True:
            try:
                self.catalog = Client.open(url, ignore_conformance=True)
                break
            except pystac_client.exceptions.APIError:
                print(f'failed opening the catalog at try {t}/{self.max_tries}')
                if t < self.max_tries:
                    t += 1
                else:
                    raise
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
    
    @staticmethod
    def _get_proc_time(scene):
        with open(os.path.join(scene, 'manifest.safe'), 'rb') as f:
            tree = etree.fromstring(f.read())
        proc = tree.find(path='.//xmlData/safe:processing',
                         namespaces=tree.nsmap)
        start = proc.attrib['start']
        del tree, proc
        return datetime.strptime(start, '%Y-%m-%dT%H:%M:%S.%f')
    
    def _filter_duplicates(self, scenes):
        tmp = sorted(scenes)
        pattern = '([0-9A-Z_]{16})_([0-9T]{15})_([0-9T]{15})'
        id_keep = []
        i = 0
        while i < len(tmp):
            group = [tmp[i]]
            match1 = re.search(pattern, os.path.basename(tmp[i])).groups()
            j = i + 1
            while j < len(tmp):
                match2 = re.search(pattern, os.path.basename(tmp[j])).groups()
                if match1 == match2:
                    group.append(tmp[j])
                    j += 1
                else:
                    break
            if len(group) > 1:
                tproc = [self._get_proc_time(x) for x in group]
                id_keep.append(i + tproc.index(max(tproc)))
            else:
                id_keep.append(i)
            i = j
        return [tmp[i] for i in id_keep]
    
    def close(self):
        del self.catalog
    
    def select(self, sensor=None, product=None, acquisition_mode=None,
               mindate=None, maxdate=None, vectorobject=None,
               date_strict=True, check_exist=True):
        """
        Select scenes from the catalog.

        Parameters
        ----------
        sensor: str or list[str] or None
            S1A or S1B
        product: str or list[str] or None
            GRD or SLC
        acquisition_mode: str or list[str] or None
            IW, EW or SM
        mindate: str or datetime.datetime or None
            the minimum acquisition date
        maxdate: str or datetime.datetime or None
            the maximum acquisition date
        vectorobject: spatialist.vector.Vector or None
            a geometry with which the scenes need to overlap
        date_strict: bool
            treat dates as strict limits or also allow flexible limits to incorporate scenes
            whose acquisition period overlaps with the defined limit?
            
            - strict: start >= mindate & stop <= maxdate
            - not strict: stop >= mindate & start <= maxdate
        check_exist: bool
            check whether found files exist locally?

        Returns
        -------
        list[str]
            the locations of the scene directories with suffix .SAFE
        
        See Also
        --------
        pystac_client.Client.search
        """
        pars = locals()
        del pars['date_strict']
        del pars['check_exist']
        del pars['self']
        
        lookup = {'product': 'sar:product_type',
                  'acquisition_mode': 'sar:instrument_mode',
                  'mindate': 'start_datetime',
                  'maxdate': 'end_datetime',
                  'sensor': 'platform'}
        lookup_platform = {'S1A': 'sentinel-1a',
                           'S1B': 'sentinel-1b'}
        
        flt = {'op': 'and', 'args': []}
        
        for key in pars.keys():
            val = pars[key]
            if val is None:
                continue
            if key == 'mindate':
                if isinstance(val, datetime):
                    val = datetime.strftime(val, '%Y%m%dT%H%M%S')
                if date_strict:
                    arg = {'op': '>=', 'args': [{'property': 'start_datetime'}, val]}
                else:
                    arg = {'op': '>=', 'args': [{'property': 'end_datetime'}, val]}
            elif key == 'maxdate':
                if isinstance(val, datetime):
                    val = datetime.strftime(val, '%Y%m%dT%H%M%S')
                if date_strict:
                    arg = {'op': '<=', 'args': [{'property': 'end_datetime'}, val]}
                else:
                    arg = {'op': '<=', 'args': [{'property': 'start_datetime'}, val]}
            elif key == 'vectorobject':
                if isinstance(val, Vector):
                    with val.clone() as vec:
                        vec.reproject(4326)
                        ext = vec.extent
                        arg = {'op': 's_intersects',
                               'args': [{'property': 'geometry'},
                                        {'type': 'Polygon',
                                         'coordinates': [[[ext['xmin'], ext['ymin']],
                                                          [ext['xmin'], ext['ymax']],
                                                          [ext['xmax'], ext['ymax']],
                                                          [ext['xmax'], ext['ymin']],
                                                          [ext['xmin'], ext['ymin']]]]}],
                               }
                else:
                    raise TypeError('argument vectorobject must be of type spatialist.vector.Vector')
            else:
                args = []
                if isinstance(val, str):
                    val = [val]
                for v in val:
                    if key == 'sensor':
                        value = lookup_platform[v]
                    else:
                        value = v
                    a = {'op': '=', 'args': [{'property': lookup[key]}, value]}
                    args.append(a)
                if len(args) == 1:
                    arg = args[0]
                else:
                    arg = {'op': 'or', 'args': args}
            flt['args'].append(arg)
        t = 1
        while True:
            try:
                result = self.catalog.search(collections=self.collections,
                                             filter=flt, max_items=None)
                result = list(result.items())
                break
            except pystac_client.exceptions.APIError:
                print(f'failed searching the catalog at try {t}/{self.max_tries}')
                if t < self.max_tries:
                    t += 1
                else:
                    raise
        out = []
        for item in result:
            assets = item.assets
            ref = assets[list(assets.keys())[0]]
            href = ref.href
            path = href[:re.search(r'\.SAFE', href).end()]
            if check_exist:
                if not Path(path).exists():
                    raise RuntimeError('scene does not exist locally:', path)
            out.append(path)
        out = self._filter_duplicates(out)
        return out
