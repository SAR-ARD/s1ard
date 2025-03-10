import re
import itertools
from lxml import html
from spatialist.vector import Vector, wkt2vector, bbox
from spatialist.auxil import utm_autodetect
from s1ard.ancillary import get_max_ext, buffer_min_overlap, get_kml, combine_polygons
from osgeo import ogr


def tile_from_aoi(vector, epsg=None, strict=True, return_geometries=False, tilenames=None):
    """
    Return a list of MGRS tile IDs or vector objects overlapping one or multiple areas of interest.
    
    Parameters
    -------
    vector: spatialist.vector.Vector or list[spatialist.vector.Vector]
        The vector object(s) to read. CRS must be EPSG:4236.
    epsg: int or list[int] or None
        Define which EPSG code(s) are allowed for the tile selection.
        If None, all tile IDs are returned regardless of projection.
    strict: bool
        Strictly only return the names/geometries of the overlapping tiles in the target projection
        or also allow reprojection of neighbouring tiles?
        In the latter case a tile name takes the form <tile ID>_<EPSG code>, e.g. `33TUN_32632`.
        Only applies if argument `epsg` is of type `int` or a list with one element.
    return_geometries: bool
        return a list of :class:`spatialist.vector.Vector` geometry objects (or just the tile names)?
    tilenames: list[str] or None
        an optional list of MGRS tile names to limit the selection
    
    Returns
    -------
    tiles: list[str or spatialist.vector.Vector]
        A list of unique MGRS tile IDs or :class:`spatialist.vector.Vector`
        objects with an attribute `mgrs` containing the tile ID.
    """
    if isinstance(epsg, int):
        epsg = [epsg]
    sortkey = None
    if return_geometries:
        sortkey = lambda x: x.mgrs
    kml = get_kml()
    with Vector(kml, driver='KML') as vec_kml:
        tiles = []
        with combine_polygons(vector, multipolygon=True) as vec_aoi:
            feature = vec_aoi.getFeatureByIndex(0)
            geom = feature.GetGeometryRef()
            vec_kml.layer.SetSpatialFilter(geom)
            feature = geom = None
            if tilenames is not None:
                values = ", ".join([f"'{x}'" for x in tilenames])
                sql_where = f"Name IN ({values})"
                vec_kml.layer.SetAttributeFilter(sql_where)
            for tile in vec_kml.layer:
                tilename = tile.GetField('Name')
                attrib = description2dict(tile.GetField('Description'))
                epsg_target = None
                if epsg is not None and attrib['EPSG'] not in epsg:
                    if len(epsg) == 1 and not strict:
                        epsg_target = int(epsg[0])
                        tilename += '_{}'.format(epsg_target)
                    else:
                        continue
                if return_geometries:
                    wkt = multipolygon2polygon(attrib['UTM_WKT'])
                    geom = wkt2vector_regrid(wkt=wkt,
                                             epsg_in=attrib['EPSG'],
                                             epsg_out=epsg_target)
                    geom.mgrs = tilename
                    tiles.append(geom)
                else:
                    tiles.append(tilename)
            tile = None
            geom = None
            layer = None
    return sorted(tiles, key=sortkey)


def aoi_from_tile(tile):
    """
    Extract one or multiple MGRS tiles from the global Sentinel-2 tiling grid and return it as a :class:`~spatialist.vector.Vector`
    object.
    
    Parameters
    ----------
    tile: str or list[str]
        The MGRS tile ID(s) that should be extracted and returned as a vector object.
        Can also be expressed as <tile ID>_<EPSG code> (e.g. `33TUN_32632`). In this case the geometry
        of the tile is reprojected to the target EPSG code, its corner coordinates rounded to multiples
        of 10, and a new :class:`~spatialist.vector.Vector` object created.
    
    Returns
    -------
    spatialist.vector.Vector or list[spatialist.vector.Vector]
        either a single object or a list depending on `tile`
    """
    kml = get_kml()
    if isinstance(tile, str):
        tile = [tile]
    
    tilenames = []
    epsg_codes = []
    pattern = '([A-Z0-9]{5})_?([0-9]+)?'
    for i in tile:
        tilename, epsg = re.search(pattern, i).groups()
        tilenames.append(tilename)
        epsg_codes.append(epsg)
    
    values = ", ".join([f"'{x}'" for x in tilenames])
    sql_where = f"Name IN ({values})"
    out = []
    with Vector(kml, driver='KML') as vec:
        vec.layer.SetAttributeFilter(sql_where)
        for i, feat in enumerate(vec.layer):
            attrib = description2dict(feat.GetField('Description'))
            wkt = multipolygon2polygon(attrib['UTM_WKT'])
            epsg_target = epsg_codes[i]
            geom = wkt2vector_regrid(wkt=wkt,
                                     epsg_in=attrib['EPSG'],
                                     epsg_out=epsg_target)
            out.append(geom)
        vec.vector.ReleaseResultSet(result_layer)
    if len(out) == 1:
        return out[0]
    else:
        return out


def description2dict(description):
    """
    Convert the HTML description field of the MGRS tile KML file to a dictionary.
    
    Parameters
    ----------
    description: str
        The plain text of the `Description` field
    
    Returns
    -------
    attrib: dict
        A dictionary with keys 'TILE_ID', 'EPSG', 'MGRS_REF', 'UTM_WKT' and 'LL_WKT'.
        The value of field 'EPSG' is of type integer, all others are strings.
    """
    attrib = html.fromstring(description)
    attrib = [x for x in attrib.xpath('//tr/td//text()') if x != ' ']
    attrib = dict(zip(attrib[0::2], attrib[1::2]))
    attrib['EPSG'] = int(attrib['EPSG'])
    return attrib


def aoi_from_scene(scene, multi=True, percent=1):
    """
    Get processing AOIs for a SAR scene. The MGRS grid requires a SAR
    scene to be geocoded to multiple UTM zones depending on the overlapping
    MGRS tiles and their projection. This function returns the following
    for each UTM zone group:
    
    - the extent in WGS84 coordinates (key `extent`)
    - the EPSG code of the UTM zone (key `epsg`)
    - the Easting coordinate for pixel alignment (key `align_x`)
    - the Northing coordinate for pixel alignment (key `align_y`)
    
    A minimum overlap of the AOIs with the SAR scene is ensured by buffering
    the AOIs if necessary. The minimum overlap can be controlled with
    parameter `percent`.
    
    Parameters
    ----------
    scene: pyroSAR.drivers.ID
        the SAR scene object
    multi: bool
        split into multiple AOIs per overlapping UTM zone or just one AOI
        covering the whole scene. In the latter case the best matching UTM
        zone is auto-detected
        (using function :func:`spatialist.auxil.utm_autodetect`).
    percent: int or float
        the minimum overlap in percent of each AOI with the SAR scene.
        See function :func:`s1ard.ancillary.buffer_min_overlap`.

    Returns
    -------
    list[dict]
        a list of dictionaries with keys `extent`, `epsg`, `align_x`,
        `align_y`
    """
    out = []
    if multi:
        # extract all overlapping tiles
        with scene.geometry() as geom:
            tiles = tile_from_aoi(vector=geom, return_geometries=True)
        
        # group tiles by UTM zone
        def fn(x):
            return x.getProjection(type='epsg')
        
        for zone, group in itertools.groupby(tiles, lambda x: fn(x)):
            geometries = list(group)
            # get UTM EPSG code
            epsg = geometries[0].getProjection(type='epsg')
            # get maximum extent of tile group
            ext = get_max_ext(geometries=geometries)
            # determine corner coordinate for alignment
            align_x = ext['xmin']
            align_y = ext['ymax']
            del geometries
            # convert extent to EPSG:4326
            with bbox(coordinates=ext, crs=epsg) as geom:
                geom.reproject(projection=4326)
                ext = geom.extent
            # ensure a minimum overlap between AOI and pre-processed scene
            with bbox(ext, 4326) as geom1:
                with scene.geometry() as geom2:
                    with buffer_min_overlap(geom1=geom1, geom2=geom2,
                                            percent=percent) as buffered:
                        ext = buffered.extent
            out.append({'extent': ext, 'epsg': epsg,
                        'align_x': align_x, 'align_y': align_y})
    else:
        with scene.bbox() as geom:
            ext = geom.extent
            # auto-detect UTM zone
            epsg = utm_autodetect(geom, 'epsg')
            # get all tiles, reprojected to the target UTM zone if necessary
            tiles = tile_from_aoi(vector=geom, epsg=epsg,
                                  return_geometries=True, strict=False)
        # determine corner coordinate for alignment
        ext_utm = tiles[0].extent
        align_x = ext_utm['xmin']
        align_y = ext_utm['ymax']
        del tiles
        out.append({'extent': ext, 'epsg': epsg,
                    'align_x': align_x, 'align_y': align_y})
    return out


def multipolygon2polygon(wkt):
    """
    Convert a MultiPolygon WKT with one Polygon to a simple Polygon WKT.
    The Sentinel-2 KML grid file stores all geometries as MultiPolygons.
    This function simply converts the geometries to simple Polygons.
    Not all geometries in the KML file have been checked. In case there are
    ever multiple Polygons in one MultiPolygon, an RuntimeError is raised.
    All other geometries are returned as is.
    
    Parameters
    ----------
    wkt: str
        A geometry WKT representation.

    Returns
    -------
    str
        the output WKT geometry. Either the original geometry or
        a Polygon extracted from a single-Polygon MultiPolygon.
    """
    geom1 = ogr.CreateGeometryFromWkt(wkt)
    poly_type = geom1.GetGeometryType()
    poly_num = geom1.GetGeometryCount()
    if poly_type == ogr.wkbMultiPolygon:
        if poly_num == 1:
            geom2 = geom1.GetGeometryRef(0).Clone()
            wkt = geom2.ExportToWkt()
            geom2 = None
        else:
            raise RuntimeError('got a MultiPolygon with multiple Polygons')
    geom1 = None
    return wkt


def wkt2vector_regrid(wkt, epsg_in, epsg_out=None):
    """
    Convert a WKT geometry to a :class:`spatialist.vector.Vector` object and
    optionally reproject and regrid it.

    Parameters
    ----------
    wkt: str
        the WKT string
    epsg_in: int
        the EPSG code for the CRS of `wkt`
    epsg_out: int or None
        and optional target CRS to reproject the geometry

    Returns
    -------
    spatialist.vector.Vector
        the geometry object
    
    See Also
    --------
    spatialist.vector.wkt2vector
    """
    if epsg_out is None:
        return wkt2vector(wkt, epsg_in)
    else:
        with wkt2vector(wkt, epsg_in) as tmp:
            tmp.reproject(epsg_out)
            ext = tmp.extent
            for k, v in ext.items():
                ext[k] = round(v / 10) * 10
        return bbox(ext, crs=epsg_out)
