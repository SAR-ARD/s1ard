S1-NRB Production
=================

The following describes the current workflow for producing the Sentinel-1 Normalised Radar Backscatter product (S1-NRB), which is being developed in study 1 of the COPA project.
This is not part of the official pyroSAR documentation.
However, as pyroSAR is the foundation of the processor, its documentation is used to outline the processor details to conveniently link to all relevant functionality.


The basis of the processing chain builds the Sentinel-2 Military Grid Reference System (MGRS) tiling system.
Hence, a reference file is needed containing the respective tile information for processing S1-NRB products.
A KML file is available online that will be used in the following steps:

https://sentinel.esa.int/documents/247904/1955685/S2A_OPER_GIP_TILPAR_MPC__20151209T095117_V20150622T000000_21000101T000000_B00.kml

This file contains all relevant information about individual tiles, in particular the EPSG code of the respective UTM zone and the geometry of the tile in UTM coordinates.
The code snippet below demonstrates the tile reading mechanism (using class :class:`spatialist.vector.Vector` and function :func:`spatialist.vector.wkt2vector`):

.. code-block:: python

    from lxml import html
    from spatialist.vector import Vector, wkt2vector

    def extract_tile(kml, tile):
        with Vector(kml, driver='KML') as vec:
            feat = vec.getFeatureByAttribute('Name', tile)
            attrib = html.fromstring(feat.GetField('Description')
            attrib = [x for x in attrib.xpath('//tr/td//text()') if x != ' ']
            attrib = dict(zip(attrib[0::2], attrib[1::2]))
            feat = None
        return wkt2vector(attrib['UTM_WKT'], int(attrib['EPSG']))

The S1 images are managed in a local SQLite database to select scenes for processing (see pyroSAR's section on `Database Handling`_).

After loading an MGRS tile as an :class:`spatialist.vector.Vector` object and selecting all relevant overlapping scenes
from the database, processing can commence.

The central function for processing backscatter data with SNAP is :func:`pyroSAR.snap.util.geocode`. It will perform all necessary steps to
generate radiometrically terrain corrected gamma naught backscatter plus all relevant additional datasets like
local incident angle and local contribution area (see argument ``export_extra``).

The code below presents an incomplete call to :func:`pyroSAR.snap.util.geocode` where several variables have been set implicitly.
``infile`` can either be  a single Sentinel-1 scene or multiple, which will then be mosaiced in radar geometry prior to geocoding.
``vec`` is the :class:`spatialist.vector.Vector` object
created from the S2 KML file with corner coordinate (``xmax``, ``ymin``). The resulting image tiles are aligned to this corner coordinate.

.. code-block:: python

    from pyroSAR.snap import geocode

    epsg = 32632  # just for demonstration; will be read from KML file
    spacing = 10
    dem = 'Copernicus 30m Global DEM'

    geocode(infile=infile, outdir=outdir, tmpdir=tmpdir,
            t_srs=epsg, shapefile=vec, tr=spacing,
            alignToStandardGrid=True,
            standardGridOriginX=xmax, standardGridOriginY=ymin,
            demName=dem, scaling='linear',
            export_extra=['localIncidenceAngle', 'incidenceAngleFromEllipsoid',
                          'scatteringArea', 'layoverShadowMask'])


The Copernicus GLO-30 DEM can easily be used for processing as it is available via SNAP auto-download. Furthermore,
Copernicus EEA-10 DEM usage has been implemented as part of the function :func:`pyroSAR.auxdata.dem_autoload`.

Many DEMs contain heights relative to a geoid such as EGM96. For SAR processing this information needs to be converted to WGS84 ellipsoid heights.
pyroSAR offers a function :func:`pyroSAR.auxdata.get_egm_lookup` to download a conversion file used by SNAP. However, SNAP itself will also automatically download this file if not found.

Alternative to the auto-download options, a custom DEM can be passed to :func:`pyroSAR.snap.util.geocode` via argument ``externalDEMFile``.
The function :func:`pyroSAR.auxdata.dem_create` can be used to directly convert between EGM96 and WGS84 heights using GDAL.
This way, the argument ``externalDEMApplyEGM`` of function :func:`pyroSAR.snap.util.geocode` can be set to ``False`` and no additional lookup file is needed.

Sentinel-1 orbit state vector files (OSV) for enhancing the orbit location accuracy are downloaded directly by pyroSAR (see :class:`pyroSAR.S1.OSV`), but can also be downloaded automatically by SNAP.
For S1-NRB processing at least Restituted Orbit files (RESORB) are needed while the more accurate Precise Orbit Ephemerides (POEORB) delivered two weeks after scene acquisition do not provide additional benefit.

The function :func:`pyroSAR.snap.util.geocode` will create a list of plain GeoTIFF files, which are slightly larger than the actual tile to ensure full tile coverage after geocoding.
These files are then subsetted to the actual tile extent, converted to Cloud Optimized GeoTIFFs (COG), and renamed to the S1-NRB naming scheme.
The function :func:`spatialist.auxil.gdalwarp` is used for this task, which is a simple wrapper around the gdalwarp utility of GDAL.
The following is another incomplete code example highlighting the general procedure of converting the individual images.
The ``outfile`` name is generated from information of the source images, the MGRS tile ID and the name of the respective file as written by :func:`pyroSAR.snap.util.geocode`.

.. code-block:: python

    from spatialist import gdalwarp, Raster
    from osgeo import gdal

    write_options = ['BLOCKSIZE=512',
                     'COMPRESS=LERC_ZSTD',
                     'MAX_Z_ERROR=0.001']

    with Raster(infiles, list_separate=False) as ras:
        source = ras.filename

    gdalwarp(src=source, dst=outfile,
             options={'format': 'COG',
                      'outputBounds': [xmin, ymin, xmax, ymax],
                      'creationOptions': write_options})

After all COG files have been created, GDAL VRT files are written for log scaling and sigma naught RTC backscatter computation.
The code below demonstrates the generation of a VRT file using :func:`spatialist.auxil.gdalbuildvrt` followed by an XML
modification to insert the pixel function (a way to achieve this with GDAL's gdalbuildvrt functionality has not yet been found).

.. code-block:: python

    from lxml import etree
    from spatialist import gdalbuildvrt

    def create_vrt(src, dst, fun, scale=None, offset=None, options=None):
        gdalbuildvrt(src=src, dst=dst, options=options)
        tree = etree.parse(dst)
        root = tree.getroot()
        band = tree.find('VRTRasterBand')
        band.attrib['subClass'] = 'VRTDerivedRasterBand'
        pixfun = etree.SubElement(band, 'PixelFunctionType')
        pixfun.text = fun
        if scale is not None:
            sc = etree.SubElement(band, 'Scale')
            sc.text = str(scale)
        if offset is not None:
            off = etree.SubElement(band, 'Offset')
            off.text = str(offset)
        etree.indent(root)
        tree.write(dst, pretty_print=True, xml_declaration=False, encoding='utf-8')

In a last step the OGC XML and STAC JSON metadata files will be written for the S1-NRB product.

.. _Database Handling: https://pyrosar.readthedocs.io/en/latest/general/processing.html#database-handling
