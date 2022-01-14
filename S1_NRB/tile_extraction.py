import os
from lxml import html
from spatialist import ogr2ogr
from spatialist.vector import Vector, wkt2vector


def tiles_from_aoi(filename, kml):
    """
    Return a list of unique MGRS tile IDs that overlap with an area of interest (AOI) provided as a vector file.
    
    Parameters
    -------
    filename: str
        The vector file to read. The following file extensions are auto-detected:
            .geojson (GeoJSON)
            .gpkg (GPKG)
            .kml (KML)
            .shp (ESRI Shapefile)
    kml: str
        Path to the Sentinel-2 tiling grid kml file provided by ESA, which can be retrieved from:
        https://sentinels.copernicus.eu/web/sentinel/missions/sentinel-2/data-products
    
    Returns
    -------
    tiles: list[str]
        A list of unique UTM tile IDs.
    """
    out = os.path.join(os.path.dirname(kml), 'tmp.gpkg')
    
    with Vector(filename) as aoi:
        if aoi.getProjection('epsg') != 4326:
            aoi.reproject(4326)
        ext = aoi.extent
        spat = (ext['xmin'], ext['ymin'], ext['xmax'], ext['ymax'])
        ogr2ogr(src=kml, dst=out, options={'format': 'GPKG', 'layers': ['Features'], 'spatFilter': spat})
    
    with Vector(out, driver='GPKG') as vec:
        tiles = vec.getUniqueAttributes('Name')
    
    os.remove(out)
    return tiles


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
    spatialist.vector.Vector object
    """
    with Vector(kml, driver='KML') as vec:
        feat = vec.getFeatureByAttribute('Name', tile)
        attrib = html.fromstring(feat.GetFieldAsString(1))
        attrib = [x for x in attrib.xpath('//tr/td//text()') if x != ' ']
        attrib = dict(zip(attrib[0::2], attrib[1::2]))
        feat = None
    return wkt2vector(attrib['UTM_WKT'], int(attrib['EPSG']))


def main(config, tr):
    """
    
    Parameters
    ----------
    config: dict
        Dictionary of the parsed config parameters for the current process.
    tr: int
        The target pixel spacing in meters, which is passed to `pyroSAR.snap.util.geocode`.
    
    Returns
    -------
    geo_dict: dict
        Dictionary containing geospatial information for each unique MGRS tile ID that will be processed.
    align_dict: dict
        Dictionary containing 'xmax' and 'ymin' coordinates to use for the `standardGridOriginX` and
        `standardGridOriginY` parameters of `pyroSAR.snap.util.geocode`
    """
    
    if config['aoi_tiles'] is not None:
        tiles = config['aoi_tiles']
    elif config['aoi_tiles'] is None and config['aoi_geometry'] is not None:
        tiles = tiles_from_aoi(filename=config['aoi_geometry'], kml=config['kml_file'])
    elif config['aoi_tiles'] is None and config['aoi_geometry'] is None:
        raise RuntimeError("Either 'aoi_tiles' or 'aoi_geometry' need to be provided!")
    
    geo_dict = {}
    for tile in tiles:
        with extract_tile(kml=config['kml_file'], tile=tile) as vec:
            ext = vec.extent
            epsg = vec.getProjection('epsg')
            xmax = ext['xmax'] - tr / 2
            ymin = ext['ymin'] + tr / 2
            geo_dict[tile] = {'ext': ext,
                              'epsg': epsg,
                              'xmax': xmax,
                              'ymin': ymin}
    
    align_dict = {'xmax': max([geo_dict[tile]['xmax'] for tile in list(geo_dict.keys())]),
                  'ymin': min([geo_dict[tile]['ymin'] for tile in list(geo_dict.keys())])}
    
    return geo_dict, align_dict
