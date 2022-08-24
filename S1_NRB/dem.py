import os
from getpass import getpass
from pyroSAR.auxdata import dem_autoload, dem_create
import S1_NRB.tile_extraction as tile_ex
from S1_NRB.ancillary import generate_unique_id
from spatialist import Raster, bbox


def prepare(geometries, dem_type, spacing, dem_dir, wbm_dir,
            kml_file, threads, epsg=None, username=None, password=None):
    """
    Downloads DEM tiles and restructures them into the MGRS tiling scheme including re-projection
    and vertical datum conversion.

    Parameters
    ----------
    geometries: list[spatialist.vector.Vector]
        A list of geometries for which to prepare the DEM tiles.
    dem_type: str
        The DEM type.
    spacing: int
        The target pixel spacing.
    dem_dir: str
        The DEM target directory.
    wbm_dir: str
        The WBM target directory.
    kml_file: str
        The KML file containing the MGRS tile geometries.
    threads: int
        The number of threads to pass to :func:`pyroSAR.auxdata.dem_create`.
    epsg: int, optional
        The coordinate reference system as an EPSG code.
    username: str or None
        The username for accessing the DEM tiles. If None and authentication is required
        for the selected DEM type, the environment variable 'DEM_USER' is read.
        If this is not set, the user is prompted interactively to provide credentials.
    password: str or None
        The password for accessing the DEM tiles.
        If None: same behavior as for username but with env. variable 'DEM_PASS'.
    """
    if dem_type == 'GETASSE30':
        geoid_convert = False
    else:
        geoid_convert = True
    geoid = 'EGM2008'
    
    buffer = 1.5  # degrees to ensure full coverage of all overlapping MGRS tiles
    tr = spacing
    wbm_dems = ['Copernicus 10m EEA DEM',
                'Copernicus 30m Global DEM II']
    dems_auth = wbm_dems
    wbm_dir = os.path.join(wbm_dir, dem_type)
    dem_dir = os.path.join(dem_dir, dem_type)
    
    if username is None:
        username = os.getenv('DEM_USER')
    if password is None:
        password = os.getenv('DEM_PASS')
    
    for i, geometry in enumerate(geometries):
        print('###### [    DEM] processing geometry {0} of {1}'.format(i + 1, len(geometries)))
        ###############################################
        tiles = tile_ex.tiles_from_aoi(vectorobject=geometry, kml=kml_file,
                                       epsg=epsg, strict=False)
        dem_names = [os.path.join(dem_dir, '{}_DEM.tif'.format(tile)) for tile in tiles]
        dem_target = {tile: name for tile, name in zip(tiles, dem_names)
                      if not os.path.isfile(name)}
        wbm_names = [os.path.join(wbm_dir, '{}_WBM.tif'.format(tile)) for tile in tiles]
        wbm_target = {tile: name for tile, name in zip(tiles, wbm_names)
                      if not os.path.isfile(name)}
        
        if len(dem_target.keys()) == 0 and len(wbm_target.keys()) == 0:
            continue
        ###############################################
        extent = geometry.extent
        ext_id = generate_unique_id(encoded_str=str(extent).encode())
        
        fname_wbm_tmp = os.path.join(wbm_dir, 'mosaic_{}.vrt'.format(ext_id))
        fname_dem_tmp = os.path.join(dem_dir, 'mosaic_{}.vrt'.format(ext_id))
        
        if not os.path.isfile(fname_wbm_tmp) or not os.path.isfile(fname_dem_tmp):
            if dem_type in dems_auth:
                if username is None:
                    username = input('Please enter your DEM access username:')
                if password is None:
                    password = getpass('Please enter your DEM access password:')
        
        print('### downloading DEM tiles')
        if dem_type in wbm_dems:
            os.makedirs(wbm_dir, exist_ok=True)
            if not os.path.isfile(fname_wbm_tmp):
                dem_autoload([geometry], demType=dem_type,
                             vrt=fname_wbm_tmp, buffer=buffer, product='wbm',
                             username=username, password=password,
                             nodata=1, hide_nodata=True)
        
        os.makedirs(dem_dir, exist_ok=True)
        if not os.path.isfile(fname_dem_tmp):
            dem_autoload([geometry], demType=dem_type,
                         vrt=fname_dem_tmp, buffer=buffer, product='dem',
                         username=username, password=password,
                         dst_nodata=0, hide_nodata=True)
        ###############################################
        if len(dem_target.keys()) > 0:
            print('### creating DEM MGRS tiles: \n{tiles}'.format(tiles=list(dem_target.keys())))
        for tilename, filename in dem_target.items():
            with tile_ex.extract_tile(kml_file, tilename) as tile:
                epsg_tile = tile.getProjection('epsg')
                ext = tile.extent
                bounds = [ext['xmin'], ext['ymin'],
                          ext['xmax'], ext['ymax']]
                dem_create(src=fname_dem_tmp, dst=filename,
                           t_srs=epsg_tile, tr=(tr, tr), pbar=True,
                           geoid_convert=geoid_convert, geoid=geoid,
                           outputBounds=bounds, threads=threads,
                           nodata=-32767)
        if os.path.isfile(fname_wbm_tmp):
            if len(wbm_target.keys()) > 0:
                print('### creating WBM MGRS tiles: \n{tiles}'.format(tiles=list(wbm_target.keys())))
            for tilename, filename in wbm_target.items():
                with tile_ex.extract_tile(kml_file, tilename) as tile:
                    epsg_tile = tile.getProjection('epsg')
                    ext = tile.extent
                    bounds = [ext['xmin'], ext['ymin'],
                              ext['xmax'], ext['ymax']]
                    dem_create(src=fname_wbm_tmp, dst=filename,
                               t_srs=epsg_tile, tr=(tr, tr),
                               resampling_method='mode', pbar=True,
                               outputBounds=bounds, threads=threads,
                               nodata='None')
        print('=' * 40)


def mosaic(geometry, dem_type, outname, epsg, kml_file, dem_dir):
    """
    Create a new mosaic GeoTIFF file from MGRS-tiled DEMs as created by :func:`S1_NRB.dem.prepare`.
    
    Parameters
    ----------
    geometry: spatialist.vector.Vector
        The geometry to be covered by the mosaic.
    dem_type: str
        The DEM type.
    outname: str
        The name of the mosaic.
    epsg: int
        The coordinate reference system as an EPSG code.
    kml_file: str
        The KML file containing the MGRS tile geometries.
    dem_dir: str
        The directory containing the DEM MGRS tiles.
    """
    dem_buffer = 200  # meters
    if not os.path.isfile(outname):
        print('### creating scene-specific DEM mosaic:', outname)
        with geometry.clone() as footprint:
            footprint.reproject(epsg)
            extent = footprint.extent
            extent['xmin'] -= dem_buffer
            extent['ymin'] -= dem_buffer
            extent['xmax'] += dem_buffer
            extent['ymax'] += dem_buffer
            with bbox(extent, epsg) as dem_box:
                tiles = tile_ex.tiles_from_aoi(vectorobject=geometry, kml=kml_file,
                                               epsg=epsg, strict=False)
                dem_names = [os.path.join(dem_dir, dem_type, '{}_DEM.tif'.format(tile)) for tile in tiles]
                with Raster(dem_names, list_separate=False)[dem_box] as dem_mosaic:
                    dem_mosaic.write(outname, format='GTiff')
