import numpy as np
from osgeo import gdal, osr
from spatialist.ancillary import finder


def extract(src, dst, variable):
    """
    Extract an OCN product's image variable and write it to a new GeoTIFF file.
    Coordinates are extracted from the corresponding latitude and longitude
    image variables and the corner coordinates written as ground control
    points (GCPs) to the output file.

    Parameters
    ----------
    src: str
        path to OCN product SAFE folder
    dst: str
        the name of the GeoTIFF file to write
    variable: str
        name of the layer to extract from the OCN product, e.g. `owiNrcsCmod`

    Returns
    -------

    """
    ocn_target = finder(target=src, matchlist=['*.nc'])[0]
    ras_ocn = gdal.Open(f'NETCDF:{ocn_target}:{variable}')
    arr = ras_ocn.ReadAsArray()
    lines, samples = arr.shape
    
    driver = gdal.GetDriverByName('GTiff')
    ras_out = driver.CreateCopy(dst, ras_ocn)
    
    ras_lat = gdal.Open(f'NETCDF:{ocn_target}:{variable[:3]}Lat')
    ras_lon = gdal.Open(f'NETCDF:{ocn_target}:{variable[:3]}Lon')
    arr_lat = ras_lat.ReadAsArray()
    arr_lon = ras_lon.ReadAsArray()
    ras_lat = ras_lon = None
    xres = (np.max(arr_lon) - np.min(arr_lon)) / samples
    yres = (np.max(arr_lat) - np.min(arr_lat)) / lines
    # coordinates are at pixel center
    ulx = float(arr_lon[0, 0]) - xres / 2
    uly = float(arr_lat[0, 0]) + yres / 2
    urx = float(arr_lon[0, -1]) - xres / 2
    ury = float(arr_lat[0, -1]) + yres / 2
    llx = float(arr_lon[-1, 0]) - xres / 2
    lly = float(arr_lat[-1, 0]) + yres / 2
    lrx = float(arr_lon[-1, -1]) - xres / 2
    lry = float(arr_lat[-1, -1]) + yres / 2
    gcps = [gdal.GCP(ulx, uly, 0, 0, 0),
            gdal.GCP(urx, ury, 0, samples, 0),
            gdal.GCP(lrx, lry, 0, samples, lines),
            gdal.GCP(llx, lly, 0, 0, lines)]
    crs = osr.SpatialReference()
    crs.SetFromUserInput('EPSG:4326')
    ras_out.SetGCPs(gcps, crs)
    outband = ras_out.GetRasterBand(1)
    outband.SetNoDataValue(-999)
    outband.WriteArray(arr, 0, 0)
    del arr, arr_lat, arr_lon
    outband = None
    ras_out = None
    ras_ocn = None


def gapfill(src, dst, md, si):
    """
    Fill gaps of an image file using GDAL.

    Parameters
    ----------
    src: str
        the source image file
    dst: str
        the destination image file with gaps filled
    md: int
        the interpolation maximum distance
    si: int
        the number of smoothing iterations

    Returns
    -------

    See Also
    --------
    osgeo.gdal.FillNodata
    """
    ras = gdal.Open(src)
    band = ras.GetRasterBand(1)
    
    driver = gdal.GetDriverByName('GTiff')
    out = driver.CreateCopy(dst, ras)
    out_band = out.GetRasterBand(1)
    arr = band.ReadAsArray()
    out_band.WriteArray(arr)
    
    mask = band.GetMaskBand()
    result = gdal.FillNodata(out_band, mask, md, si)
    
    out_band = None
    out = None
    band = None
    ras = None
