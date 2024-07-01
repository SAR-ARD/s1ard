import os
import re
import tempfile
import itertools
from getpass import getpass
from pyroSAR.auxdata import dem_autoload, dem_create
import s1ard.tile_extraction as tile_ex
from s1ard.ancillary import generate_unique_id, get_max_ext, vrt_add_overviews
from spatialist import Raster, bbox
import logging

log = logging.getLogger('s1ard')


def prepare(vector, dem_type, dem_dir, wbm_dir, dem_strict=True,
            tilenames=None, threads=None, username=None, password=None):
    """
    Downloads DEM and WBM tiles and restructures them into the MGRS tiling
    scheme including re-projection and vertical datum conversion.

    Parameters
    ----------
    vector: spatialist.vector.Vector
        The vector object for which to prepare the DEM and WBM tiles.
        CRS must be EPSG:4236.
    dem_type: str
        The DEM type.
    dem_dir: str or None
        The DEM target directory. DEM preparation can be skipped if set to None.
    wbm_dir: str or None
        The WBM target directory. WBM preparation can be skipped if set to None
    dem_strict: bool
        strictly only create DEM tiles in the native CRS of the MGRS tile or
        also allow reprojection to ensure full coverage of the vector object in every CRS.
    tilenames: list[str] or None
        an optional list of MGRS tile names. Default None: process all overalapping tiles.
    threads: int or None
        The number of threads to pass to :func:`pyroSAR.auxdata.dem_create`.
        Default `None`: use the value of `GDAL_NUM_THREADS` without modification.
    username: str or None
        The username for accessing the DEM tiles. If None and authentication is required
        for the selected DEM type, the environment variable 'DEM_USER' is read.
        If this is not set, the user is prompted interactively to provide credentials.
    password: str or None
        The password for accessing the DEM tiles.
        If None: same behavior as for username but with env. variable 'DEM_PASS'.
    
    Examples
    --------
    >>> from s1ard import dem
    >>> from spatialist import bbox
    >>> ext = {'xmin': 12, 'xmax': 13, 'ymin': 50, 'ymax': 51}
    # strictly only create overlapping DEM tiles in their native CRS.
    # Will create tiles 32UQA, 32UQB, 33UUR and 33UUS.
    >>> with bbox(coordinates=ext, crs=4326) as vec:
    >>>     dem.prepare(vector=vec, dem_type='Copernicus 30m Global DEM',
    >>>                 dem_dir='DEM', wbm_dir=None, dem_strict=True,
    >>>                 threads=4)
    # Process all overlapping DEM tiles to each CRS.
    # Will additionally create tiles 32UQA_32633, 32UQB_32633, 33UUR_32632 and 33UUS_32632.
    >>> with bbox(coordinates=ext, crs=4326) as vec:
    >>>     dem.prepare(vector=vec, dem_type='Copernicus 30m Global DEM',
    >>>                 dem_dir='DEM', wbm_dir=None, dem_strict=False,
    >>>                 threads=4)
    
    See Also
    --------
    s1ard.tile_extraction.tile_from_aoi
    """
    if dem_type == 'GETASSE30':
        geoid_convert = False
    else:
        geoid_convert = True
    geoid = 'EGM2008'  # applies to all Copernicus DEM options
    
    tr = 10  # target resolution. Lower resolutions can be created virtually using VRTs.
    # additional creation options for gdalwarp
    create_options = ['COMPRESS=LERC_ZSTD', 'MAX_Z_ERROR=0']
    
    # DEM options with WBMs
    wbm_dems = ['Copernicus 10m EEA DEM',
                'Copernicus 30m Global DEM',
                'Copernicus 30m Global DEM II']
    if wbm_dir is not None and dem_type in wbm_dems:
        wbm_dir = os.path.join(wbm_dir, dem_type)
    else:
        wbm_dir = None
    if dem_dir is not None:
        dem_dir = os.path.join(dem_dir, dem_type)
    if wbm_dir is None and dem_dir is None:
        return
    # get the geometries of all tiles overlapping with the AOI
    tiles = tile_ex.tile_from_aoi(vector=vector,
                                  return_geometries=True,
                                  tilenames=tilenames)
    # group the returned tiles by CRS and process them separately
    for epsg, group in itertools.groupby(tiles, lambda x: x.getProjection('epsg')):
        vectors = list(group)
        
        # In case the DEM tiles are to be prepared as well, create a new list of tiles.
        # This new list contains all tiles covering the AOI but all re-projected to the
        # current CRS. This way, a DEM mosaic can be created from prepared tiles in any
        # UTM zone covering the AOI while fully covering it. This was needed for processing
        # full SAR scenes to different UTM zones. In the current workflow this in no longer
        # used.
        if dem_dir is not None and not dem_strict:
            vectors = tile_ex.tile_from_aoi(
                vector=vector.bbox(),
                epsg=epsg,
                strict=False,
                return_geometries=True,
                tilenames=tilenames)
        
        # Get the bounding box of the tile vector objects and use this from here on
        ext = get_max_ext(geometries=vectors, buffer=200)
        with bbox(coordinates=ext, crs=epsg) as box:
            box.reproject(4326)
            ext_4326 = box.extent
        
        # create a unique ID for the names of the VRTs that will be created.
        ext_id = generate_unique_id(encoded_str=str(ext_4326).encode())
        
        if dem_dir is not None:
            dem_names_base = ['{}_DEM.tif'.format(tile.mgrs) for tile in vectors]
            dem_names = [os.path.join(dem_dir, x) for x in dem_names_base]
            dem_target = [(tile, name) for tile, name in zip(vectors, dem_names)
                          if not os.path.isfile(name)]
            fname_dem_tmp = os.path.join(dem_dir, 'mosaic_{}.vrt'.format(ext_id))
        else:
            dem_target = dict()
            fname_dem_tmp = None
        if wbm_dir is not None:
            # exclude the reprojected tiles from the list of WBM tiles
            tiles_wbm = [x for x in vectors if not re.search('_[0-9]*', x.mgrs)]
            wbm_names_base = ['{}_WBM.tif'.format(tile.mgrs) for tile in tiles_wbm]
            wbm_names = [os.path.join(wbm_dir, x) for x in wbm_names_base]
            wbm_target = [(tile, name) for tile, name in zip(tiles_wbm, wbm_names)
                          if not os.path.isfile(name)]
            fname_wbm_tmp = os.path.join(wbm_dir, 'mosaic_{}.vrt'.format(ext_id))
        else:
            wbm_target = dict()
            fname_wbm_tmp = None
        
        # stop if no files need to be created
        if len(dem_target) == 0 and len(wbm_target) == 0:
            continue
        ###############################################
        # DEM download and VRT mosaic creation
        
        # get download authentication if either WBM or DEM VRTs will be created
        c_wbm = fname_wbm_tmp is not None and not os.path.isfile(fname_wbm_tmp)
        c_dem = fname_dem_tmp is not None and not os.path.isfile(fname_dem_tmp)
        if c_wbm or c_dem:
            username, password = authenticate(dem_type=dem_type,
                                              username=username,
                                              password=password)
        
        # download WBM tiles and combine them in a VRT mosaic
        if c_wbm:
            os.makedirs(wbm_dir, exist_ok=True)
            with bbox(coordinates=ext_4326, crs=4326) as vec:
                dem_autoload(geometries=[vec], demType=dem_type,
                             vrt=fname_wbm_tmp, product='wbm',
                             username=username, password=password,
                             crop=False)
        # download DEM tiles and combine them in a VRT mosaic
        if c_dem:
            os.makedirs(dem_dir, exist_ok=True)
            with bbox(coordinates=ext_4326, crs=4326) as vec:
                dem_autoload(geometries=[vec], demType=dem_type,
                             vrt=fname_dem_tmp, product='dem',
                             username=username, password=password,
                             crop=False)
        ###############################################
        if len(dem_target) > 0:
            tiles = [x[0].mgrs for x in dem_target]
            log.info(f'creating DEM MGRS tiles: {tiles}')
        for tile, filename in dem_target:
            ext = tile.extent
            bounds = [ext['xmin'], ext['ymin'],
                      ext['xmax'], ext['ymax']]
            dem_create(src=fname_dem_tmp, dst=filename,
                       t_srs=epsg, tr=(tr, tr), pbar=False,
                       geoid_convert=geoid_convert, geoid=geoid,
                       outputBounds=bounds, threads=threads,
                       nodata=-32767, creationOptions=create_options)
        ###############################################
        if len(wbm_target) > 0:
            tiles = [x[0].mgrs for x in wbm_target]
            log.info(f'creating WBM MGRS tiles: {tiles}')
        for tile, filename in wbm_target:
            ext = tile.extent
            bounds = [ext['xmin'], ext['ymin'],
                      ext['xmax'], ext['ymax']]
            dem_create(src=fname_wbm_tmp, dst=filename,
                       t_srs=epsg, tr=(tr, tr),
                       resampleAlg='mode', pbar=False,
                       outputBounds=bounds, threads=threads,
                       creationOptions=create_options)


def authenticate(dem_type, username=None, password=None):
    """
    Query the username and password. If None, environment variables DEM_USER and DEM_PASS are read.
    If they are also None, the user is queried interactively.
    
    Parameters
    ----------
    dem_type: str
        the DEM type. Needed for determining whether authentication is needed.
    username: str or None
        The username for accessing the DEM tiles. If None and authentication is required
        for the selected DEM type, the environment variable 'DEM_USER' is read.
        If this is not set, the user is prompted interactively to provide credentials.
    password: str or None
        The password for accessing the DEM tiles.
        If None: same behavior as for username but with env. variable 'DEM_PASS'.

    Returns
    -------
    tuple[str or None]
        the username and password
    """
    dems_auth = ['Copernicus 10m EEA DEM',
                 'Copernicus 30m Global DEM II']
    if dem_type not in dems_auth:
        return None, None
    
    if username is None:
        username = os.getenv('DEM_USER')
    if password is None:
        password = os.getenv('DEM_PASS')
    if username is None:
        username = input('Please enter your DEM access username:')
    if password is None:
        password = getpass('Please enter your DEM access password:')
    return username, password


def mosaic(geometry, dem_type, outname, epsg=None, kml_file=None,
           dem_dir=None, username=None, password=None, threads=4):
    """
    Create a new scene-specific DEM mosaic GeoTIFF file.
    Can be created from MGRS-tiled DEMs as created by :func:`s1ard.dem.prepare`
    or ad hoc using :func:`pyroSAR.auxdata.dem_autoload` and :func:`pyroSAR.auxdata.dem_create`.
    In the former case the arguments `username`, `password` and `threads` are ignored and
    all tiles found in `dem_dir` are read.
    In the latter case the arguments `epsg`, `kml_file` and `dem_dir` are ignored and the DEM is
    only mosaiced and geoid-corrected.
    
    Parameters
    ----------
    geometry: spatialist.vector.Vector
        The geometry to be covered by the mosaic.
    dem_type: str
        The DEM type.
    outname: str
        The name of the mosaic.
    epsg: int or None
        The coordinate reference system as an EPSG code.
    kml_file: str or None
        The KML file containing the MGRS tile geometries.
    dem_dir: str or None
        The directory containing the DEM MGRS tiles.
    username: str or None
        The username for accessing the DEM tiles. If None and authentication is required
        for the selected DEM type, the environment variable 'DEM_USER' is read.
        If this is not set, the user is prompted interactively to provide credentials.
    password: str or None
        The password for accessing the DEM tiles.
        If None: same behavior as for username but with env. variable 'DEM_PASS'.
    threads: int
        The number of threads to pass to :func:`pyroSAR.auxdata.dem_create`.
    """
    if not os.path.isfile(outname):
        if dem_dir is not None:
            dem_buffer = 200  # meters
            with geometry.clone() as footprint:
                footprint.reproject(epsg)
                extent = footprint.extent
                extent['xmin'] -= dem_buffer
                extent['ymin'] -= dem_buffer
                extent['xmax'] += dem_buffer
                extent['ymax'] += dem_buffer
                with bbox(extent, epsg) as dem_box:
                    tiles = tile_ex.tile_from_aoi(vector=geometry, kml=kml_file,
                                                  epsg=epsg, strict=False)
                    dem_names = [os.path.join(dem_dir, dem_type, '{}_DEM.tif'.format(tile)) for tile in tiles]
                    with Raster(dem_names, list_separate=False)[dem_box] as dem_mosaic:
                        dem_mosaic.write(outname, format='GTiff')
        else:
            username, password = authenticate(dem_type=dem_type, username=username, password=password)
            buffer = 0.01  # degrees
            if dem_type == 'GETASSE30':
                geoid_convert = False
            else:
                geoid_convert = True
            geoid = 'EGM2008'
            vrt = outname.replace('.tif', '.vrt')
            dem_autoload([geometry], demType=dem_type,
                         vrt=vrt, buffer=buffer, product='dem',
                         username=username, password=password)
            dem_create(src=vrt, dst=outname, pbar=False,
                       geoid_convert=geoid_convert, geoid=geoid,
                       threads=threads, nodata=-32767)


def to_mgrs(tile, dst, dem_type, overviews, tr, format='COG',
            create_options=None, threads=None, pbar=False):
    """
    Create an MGRS-tiled DEM file.
    
    Parameters
    ----------
    tile: str
        the MGRS tile ID
    dst: str
        the destination file name
    dem_type: str
        The DEM type.
    overviews: list[int]
        The overview levels
    tr: tuple[int or float]
        the target resolution as (x, y)
    format: str
        the output file format
    create_options: list[str] or None
        additional creation options to be passed to :func:`spatialist.auxil.gdalwarp`.
    threads: int or None
        The number of threads to pass to :func:`pyroSAR.auxdata.dem_create`.
        Default `None`: use the value of `GDAL_NUM_THREADS` without modification.
    pbar: bool

    Returns
    -------

    """
    if dem_type == 'GETASSE30':
        geoid_convert = False
    else:
        geoid_convert = True
    geoid = 'EGM2008'  # applies to all Copernicus DEM options
    with tile_ex.aoi_from_tile(tile=tile) as vec:
        ext = vec.extent
        epsg = vec.getProjection('epsg')
    bounds = [ext['xmin'], ext['ymin'], ext['xmax'], ext['ymax']]
    buffer = 200
    ext['xmin'] -= buffer
    ext['ymin'] -= buffer
    ext['xmax'] += buffer
    ext['ymax'] += buffer
    vrt = tempfile.NamedTemporaryFile(suffix='.vrt').name
    with bbox(coordinates=ext, crs=epsg) as vec:
        vec.reproject(4326)
        dem_autoload(geometries=[vec], demType=dem_type, vrt=vrt)
    vrt_add_overviews(vrt=vrt, overviews=overviews)
    dem_create(src=vrt, dst=dst, t_srs=epsg, tr=tr,
               geoid_convert=geoid_convert, geoid=geoid, pbar=pbar,
               outputBounds=bounds, threads=threads, format=format,
               creationOptions=create_options)
