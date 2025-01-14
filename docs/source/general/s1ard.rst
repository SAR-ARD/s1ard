ARD Production
==============

The following sections give a brief overview of the major components of creating a S1-NRB product.
All steps are comprised in function :func:`s1ard.processor.main`.
The pyroSAR package builds the foundation of the processor and its documentation is used to outline the processor details to conveniently link to all relevant functionality.

MGRS Gridding
-------------

The basis of the processing chain builds the Sentinel-2 Military Grid Reference System (MGRS) tiling system.
Hence, a reference file is needed containing the respective tile information for processing S1-NRB products.
A KML file is available online that will be used in the following steps:

`S2A_OPER_GIP_TILPAR_MPC__20151209T095117_V20150622T000000_21000101T000000_B00.kml <https://sentinel.esa.int/documents/247904/1955685/S2A_OPER_GIP_TILPAR_MPC__20151209T095117_V20150622T000000_21000101T000000_B00.kml>`_

This file contains all relevant information about individual tiles, in particular the EPSG code of the respective UTM zone and the geometry of the tile in UTM coordinates.
This file is automatically downloaded to `~/s1ard` by the function :func:`s1ard.ancillary.get_kml`.
The function :func:`s1ard.tile_extraction.aoi_from_tile` can be used to extract one or multiple tiles as :class:`spatialist.vector.Vector` object.

Scene Management
----------------

The S1 images are managed in a local SQLite database to select scenes for processing (see pyroSAR's section on `Database Handling`_) or are directly queried from a STAC catalog (see :class:`s1ard.archive.STACArchive`).
See documentation section :doc:`/general/search` for details.

After loading an MGRS tile as an :class:`spatialist.vector.Vector` object and selecting all relevant overlapping scenes
from the database, processing can commence.

DEM Handling
------------

s1ard offers a convenience function :func:`s1ard.dem.mosaic` for creating scene-specific DEM files from various sources.
The function is based on :func:`pyroSAR.auxdata.dem_autoload` and :func:`pyroSAR.auxdata.dem_create` and will

- download all tiles of the selected source overlapping with a defined geometry
- create a GDAL VRT virtual mosaic from the tiles including gap filling over ocean areas
- create a new GeoTIFF from the VRT including geoid-ellipsoid height conversion if necessary
  (WGS84 heights are generally required for SAR processing but provided heights might be relative to a geoid like EGM2008).

OSV Handling
------------

Sentinel-1 orbit state vector files (OSV) for enhancing the orbit location accuracy are downloaded directly by pyroSAR (see :class:`pyroSAR.S1.OSV`), but can also be downloaded automatically by SNAP.
For S1-NRB processing at least Restituted Orbit files (RESORB) are needed while the more accurate Precise Orbit Ephemerides (POEORB) delivered two weeks after scene acquisition do not provide much additional benefit.

SNAP Processing
---------------

The central function for processing backscatter data with SNAP is :func:`s1ard.snap.process`. It will perform all necessary steps to
generate radiometrically terrain corrected gamma/sigma naught backscatter plus all relevant additional datasets like
local incident angle and local contribution area (see argument ``export_extra``).
In a full processor run, the following functions are called in sequence:

- :func:`s1ard.snap.pre`: general pre-processing including

  + Orbit state vector enhancement
  + (GRD only) border noise removal
  + Calibration to beta naught (for RTC) and sigma naught (for NESZ)
  + Thermal noise removal (including generation of noise equivalent sigma zero (NESZ) noise power images)
  + (SLC only) debursting and swath merging

- :func:`s1ard.snap.mli`: creates multi-looked image files (MLIs) per polarization if the target pixel spacing is larger than the source pixel spacing.

- :func:`s1ard.snap.rtc`: radiometric terrain flattening.
  Output is backscatter in gamma naught RTC (:math:`\gamma^0_T`) and sigma naught RTC (:math:`\sigma^0_T`) as well as the scattering area (:math:`\beta^0 / \gamma^0_T`).

- :func:`s1ard.snap.gsr`: computation of the gamma-sigma ratio (:math:`\sigma^0_T / \gamma^0_T`).

- :func:`s1ard.snap.geo`: geocoding. This function may be called multiple times if the scene overlaps with multiple UTM zones.

The output is a BEAM-DIMAP product which consists of a `dim` metadata file and a `data` folder containing the individual image layers in ENVI format (extension `img`).
The function :func:`s1ard.snap.find_datasets` can be used to collect the individual images files for a scene.

Depending on the user configuration parameters ``measurement`` and ``annotation``, some modifications to the workflow above are possible:

- :func:`s1ard.snap.gsr` may be replaced by :func:`s1ard.snap.sgr` to create a sigma-gamma ratio (:math:`\gamma^0_T / \sigma^0_T`)

ARD Formatting
--------------

During SAR processing, files covering a whole scene are created. In this last step, the scene-based structure is converted to the MGRS tile structure.
If one tile overlaps with multiple scenes, these scenes are first virtually mosaiced using VRT files.
The files are then subsetted to the actual tile extent, converted to Cloud Optimized GeoTIFFs (COG), and renamed to the S1-NRB or S1-ORB naming scheme.
All steps are performed by :func:`s1ard.nrb.format`.
The actual file format conversion is done with :func:`spatialist.auxil.gdalwarp`, which is a simple wrapper around the gdalwarp utility of GDAL.
The following is an incomplete code example highlighting the general procedure of converting the individual images.
The ``outfile`` name is generated from information of the source images, the MGRS tile ID and the name of the respective file of the SAR processing step.

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

After all COG files have been created, GDAL VRT files are written for log scaling and conversion to other backscatter conventions using function :func:`s1ard.nrb.create_vrt`.
The code below demonstrates the generation of a VRT file for log-scaling using :func:`spatialist.auxil.gdalbuildvrt` followed by an XML
modification to insert the pixel function (a way to achieve this with GDAL's gdalbuildvrt functionality has not yet been found).

.. code-block:: python

    from lxml import etree
    from spatialist import gdalbuildvrt

    src = 'test.tif'
    dst = 'test_db.vrt'

    gdalbuildvrt(src=src, dst=dst)
    tree = etree.parse(dst)
    root = tree.getroot()
    band = tree.find('VRTRasterBand')
    band.attrib['subClass'] = 'VRTDerivedRasterBand'
    pixfun = etree.SubElement(band, 'PixelFunctionType')
    pixfun.text = 'dB'
    arg = etree.SubElement(band, 'PixelFunctionArguments')
    arg.attrib['fact'] = '10'
    etree.indent(root)
    tree.write(dst, pretty_print=True, xml_declaration=False, encoding='utf-8')

In a last step the OGC XML and STAC JSON metadata files will be written for the S1-NRB product.

.. _Database Handling: https://pyrosar.readthedocs.io/en/latest/general/processing.html#database-handling
