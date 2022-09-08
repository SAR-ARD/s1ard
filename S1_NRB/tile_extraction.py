import re
from lxml import html
from spatialist.vector import Vector, wkt2vector, bbox
from S1_NRB.ancillary import get_max_ext


def tiles_from_aoi(vectorobject, kml, epsg=None, strict=True):
    """
    Return a list of unique MGRS tile IDs that overlap with an area of interest (AOI) provided as a
    :class:`~spatialist.vector.Vector` object.
    
    Parameters
    -------
    vectorobject: spatialist.vector.Vector
        The vector object to read.
    kml: str
        Path to the Sentinel-2 tiling grid KML file.
    epsg: int or list[int] or None
        Define which EPSG code(s) are allowed for the tile selection.
        If None, all tile IDs are returned regardless of projection.
    strict: bool
        Strictly only return the names of the overlapping tiles in the target projection
        or also allow reprojection of neighbouring tiles?
        In the latter case a tile name takes the form <tile ID>_<EPSG code>, e.g. `33TUN_32632`.
        Only applies if argument `epsg` is of type `int` or a list with one element.
    
    Returns
    -------
    tiles: list[str]
        A list of unique MGRS tile IDs.
    
    Notes
    -----
    The global Sentinel-2 tiling grid can be retrieved from:
    https://sentinel.esa.int/documents/247904/1955685/S2A_OPER_GIP_TILPAR_MPC__20151209T095117_V20150622T000000_21000101T000000_B00.kml
    """
    if isinstance(epsg, int):
        epsg = [epsg]
    if vectorobject.getProjection('epsg') != 4326:
        raise RuntimeError('the CRS of the input vector object must be EPSG:4326')
    with Vector(kml, driver='KML') as vec:
        tilenames = []
        vectorobject.layer.ResetReading()
        for item in vectorobject.layer:
            geom = item.GetGeometryRef()
            vec.layer.SetSpatialFilter(geom)
            for tile in vec.layer:
                tilename = tile.GetField('Name')
                if tilename not in tilenames:
                    attrib = description2dict(tile.GetField('Description'))
                    if epsg is not None and attrib['EPSG'] not in epsg:
                        if len(epsg) == 1 and not strict:
                            tilename += '_{}'.format(epsg[0])
                        else:
                            continue
                    tilenames.append(tilename)
        vectorobject.layer.ResetReading()
        tile = None
        geom = None
        item = None
        return sorted(tilenames)


def aoi_from_tiles(kml, tiles):
    """
    Returns the bounding box of a list of MGRS tile IDs as a :class:`~spatialist.vector.Vector` object.
    
    Parameters
    ----------
    kml: str
        Path to the Sentinel-2 tiling grid KML file.
    tiles: list[str]
        A list of unique MGRS tile IDs.
    
    Returns
    -------
    spatialist.vector.Vector
    
    Notes
    -----
    The global Sentinel-2 tiling grid can be retrieved from:
    https://sentinel.esa.int/documents/247904/1955685/S2A_OPER_GIP_TILPAR_MPC__20151209T095117_V20150622T000000_21000101T000000_B00.kml
    """
    geometries = []
    for tile in tiles:
        geom = extract_tile(kml=kml, tile=tile)
        geom.reproject(4326)
        geometries.append(geom)
    max_ext = get_max_ext(geometries=geometries)
    return bbox(max_ext, crs=4326)


def extract_tile(kml, tile):
    """
    Extract an MGRS tile from the global Sentinel-2 tiling grid and return it as a :class:`~spatialist.vector.Vector`
    object.
    
    Parameters
    ----------
    kml: str
        Path to the Sentinel-2 tiling grid KML file.
    tile: str
        The MGRS tile ID that should be extracted and returned as a vector object.
        Can also be expressed as <tile ID>_<EPSG code> (e.g. `33TUN_32632`). In this case the geometry
        of the tile is reprojected to the target EPSG code, its corner coordinates rounded to multiples
        of 10, and a new :class:`~spatialist.vector.Vector` object created.
    
    Returns
    -------
    spatialist.vector.Vector
    
    Notes
    -----
    The global Sentinel-2 tiling grid can be retrieved from:
    https://sentinel.esa.int/documents/247904/1955685/S2A_OPER_GIP_TILPAR_MPC__20151209T095117_V20150622T000000_21000101T000000_B00.kml
    """
    tilename, epsg = re.search('([A-Z0-9]{5})_?([0-9]+)?', tile).groups()
    with Vector(kml, driver='KML') as vec:
        feat = vec.getFeatureByAttribute('Name', tilename)
        attrib = description2dict(feat.GetField('Description'))
        feat = None
    if epsg is None:
        return wkt2vector(attrib['UTM_WKT'], attrib['EPSG'])
    else:
        with wkt2vector(attrib['UTM_WKT'], attrib['EPSG']) as tmp:
            tmp.reproject(int(epsg))
            ext = tmp.extent
            for k, v in ext.items():
                ext[k] = round(v / 10) * 10
        return bbox(ext, crs=int(epsg))


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


def get_tile_dict(kml_file, spacing, aoi_geometry=None, aoi_tiles=None):
    """
    Creates a dictionary with information for each unique MGRS tile ID that is being processed (extent, epsg code) as
    well as alignment coordinates that can be passed to the `standardGridOriginX` and `standardGridOriginY` parameters
    of :func:`pyroSAR.snap.util.geocode`
    
    Parameters
    ----------
    kml_file: str
        The KML file containing the MGRS tile geometries.
    spacing: int
        The target pixel spacing in meters, which is passed to :func:`pyroSAR.snap.util.geocode`.
    aoi_geometry: str or None
        A vector geometry file name.
    aoi_tiles: list[str] or None
        a list with MGRS tile names.
    
    Returns
    -------
    tile_dict: dict
        The output dictionary containing information about each unique MGRS tile ID and alignment coordinates.
    """
    if aoi_geometry is not None:
        with Vector(aoi_geometry) as aoi:
            tiles = tiles_from_aoi(aoi, kml=kml_file)
    elif aoi_tiles is not None:
        tiles = aoi_tiles
    else:
        raise RuntimeError("either 'aoi_geometry' or 'aoi_tiles' must be defined")
    
    geo_dict = {}
    for tile in tiles:
        with extract_tile(kml=kml_file, tile=tile) as vec:
            ext = vec.extent
            epsg = vec.getProjection('epsg')
            xmax = ext['xmax'] - spacing / 2
            ymin = ext['ymin'] + spacing / 2
            geo_dict[tile] = {'ext': ext,
                              'epsg': epsg,
                              'xmax': xmax,
                              'ymin': ymin}
    
    align_dict = {'xmax': max([geo_dict[tile]['xmax'] for tile in list(geo_dict.keys())]),
                  'ymin': min([geo_dict[tile]['ymin'] for tile in list(geo_dict.keys())])}
    geo_dict['align'] = align_dict
    
    return geo_dict
