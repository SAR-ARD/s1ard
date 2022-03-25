from lxml import html
from spatialist.vector import Vector, wkt2vector


def tiles_from_aoi(vectorobject, kml, epsg=None):
    """
    Return a list of unique MGRS tile IDs that overlap with an area of interest (AOI) provided as a vector object.
    
    Parameters
    -------
    vectorobject: spatialist.vector.Vector
        The vector object to read.
    kml: str
        Path to the Sentinel-2 tiling grid kml file provided by ESA, which can be retrieved from:
        https://sentinels.copernicus.eu/web/sentinel/missions/sentinel-2/data-products
    epsg: int or list[int]
        define which EPSG code(s) are allowed for the tile selection.
    
    Returns
    -------
    tiles: list[str]
        A list of unique UTM tile IDs.
    """
    if isinstance(epsg, int):
        epsg = [epsg]
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
                        continue
                    tilenames.append(tilename)
        vectorobject.layer.ResetReading()
        tile = None
        geom = None
        item = None
        return sorted(tilenames)


def extract_tile(kml, tile):
    """
    Extract a MGRS tile from the global Sentinel-2 tiling grid and return it as a vector object.
    
    Parameters
    ----------
    kml: str
        Path to the Sentinel-2 tiling grid kml file provided by ESA, which can be retrieved from:
        https://sentinels.copernicus.eu/web/sentinel/missions/sentinel-2/data-products
    tile: str
        The MGRS tile ID that should be extracted and returned as a vector object.
    
    Returns
    -------
    spatialist.vector.Vector
    """
    with Vector(kml, driver='KML') as vec:
        feat = vec.getFeatureByAttribute('Name', tile)
        attrib = description2dict(feat.GetField('Description'))
        feat = None
    return wkt2vector(attrib['UTM_WKT'], attrib['EPSG'])


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


def get_tile_dict(config, spacing):
    """
    Creates a dictionary with information for each unique MGRS tile ID that is being processed (extent, epsg code) as
    well as alignment coordinates that can be passed to the `standardGridOriginX` and `standardGridOriginY` parameters
    of `pyroSAR.snap.util.geocode`
    
    Parameters
    ----------
    config: dict
        Dictionary of the parsed config parameters for the current process.
    spacing: int
        The target pixel spacing in meters, which is passed to `pyroSAR.snap.util.geocode`.
    
    Returns
    -------
    tile_dict: dict
        The output dictionary containing information about each unique MGRS tile ID and alignment coordinates.
    """
    try:
        with Vector(config['aoi_geometry']) as aoi:
            tiles = tiles_from_aoi(aoi, kml=config['kml_file'])
    except AttributeError:
        tiles = config['aoi_tiles']
    
    geo_dict = {}
    for tile in tiles:
        with extract_tile(kml=config['kml_file'], tile=tile) as vec:
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
