Usage
=====

This section outlines how to configure and run the processor. Configuration is most conveniently kept in a ``config.ini``
configuration file but can also be modified via the command line.

Two different types of product were intended when developing the processor, Normalised Radar Backscatter (NRB)
and Ocean Radar Backscatter (ORB). However, the processor does not strictly separate between them and products can be created that
conform to both types.

To create an NRB product as defined by the CEOS ARD specification, the following configuration would be necessary:

.. code-block:: ini

    mode = sar, nrb
    measurement = gamma
    annotation = dm, ei, em, id, lc, li, np, ratio

The generated backscatter is gamma nought RTC. Annotation layers are a data mask, ellipsoidal incident angle, elevation model,
acquisition id mask, local contributing area, local incidence angle, noise power and a gamma-sigma ratio.

For ORB, the following configuration is foreseen:

.. code-block:: ini

    mode = sar, orb
    measurement = sigma
    annotation = dm, id, ld, li, np, wm

Compared to NRB, the backscatter is now sigma nought RTC. The ellipsoidal incident angle is excluded because over ocean
it is nearly identical to the local incident angle. Furthermore, the elevation model, local contributing area and
backscatter ratio are excluded as they are not seen as necessary.
Two new annotation layers are added. A look direction angle and a wind model.

See below for further details.

Configuration
-------------
Usage of the s1ard package relies on a configuration file that needs to be set up by the user. The configuration
file follows the INI format, which uses plain text to store properties as key-value pairs. INI files can be created and
opened with any text editor. An example ``config.ini`` file for the s1ard package can be found here:

https://github.com/SAR-ARD/s1ard/blob/main/config.ini

Configuration files in INI format can have different sections. Each section begins at a section name and ends at the next
section name. The ``config.ini`` file used with the s1ard package should at least have a dedicated section for processing
related parameters. This section is by default named ``[PROCESSING]``.
Users might create several processing sections in the same configuration file with parameter values that correspond to different
processing scenarios (e.g., for different areas of interest). Note that each section must contain all necessary
configuration parameters even if only a few are varied between the sections.

The following provides an overview of the parameters the ``config.ini`` should contain and anything that should be
considered when selecting their values:

Processing Section
^^^^^^^^^^^^^^^^^^

mode
++++

Options: ``sar | nrb | orb``

This parameter determines what steps should be executed.
``sar`` will only start SAR preprocessing, whereas ``nrb`` and ``orb`` will only start ARD generation from existing SAR 
products preprocessed in ``sar``.
By defining both ``sar`` and one of the ARD modes as list, both SAR preprocessing and ARD generation can be run together:

.. code-block:: python

    mode = sar, nrb

scene
+++++

Define a single SAR scene filename instead of searching for scenes in a database.
If this parameter is set, the ``mode`` must be ``sar``.
In case of a GRD, database search is still performed to collect neighbors.

aoi_tiles & aoi_geometry
++++++++++++++++++++++++

Limit processing to a specific area of interest (AOI).

``aoi_tiles`` can be used to define the area of interest via MGRS tile IDs, which must be provided comma-separated (e.g.,
``aoi_tiles = 32TNS, 32TMT, 32TMS``). ``aoi_geometry`` defines the area of interest via a full path to a vector file
supported by :class:`spatialist.vector.Vector`. This option will automatically search for overlapping MGRS tiles and use
these for processing.
Both parameters are optional and can be set to ``None`` or left empty. ``aoi_tiles`` overrides ``aoi_geometry``.
If neither is defined, all tiles overlapping with the scene search result are processed.

mindate & maxdate
+++++++++++++++++

Search for source scenes within the defined date range.
Allowed are all string representations that can be parsed by :meth:`dateutil.parser.parse`.

date_strict
+++++++++++

Treat dates as strict limits or also allow flexible limits to incorporate scenes
whose acquisition period overlaps with the defined limit.

 + strict: ``start >= mindate & stop <= maxdate``
 + not strict: ``stop >= mindate & start <= maxdate``

sensor
++++++

Options: ``S1A | S1B``

The Sentinel-1 sensor/platform.

acq_mode
++++++++

Options: ``IW | EW | SM``

The acquisition mode of the source scenes that should be processed.

product
+++++++

Options: ``GRD | SLC``

The product of the source scenes that should be processed.

datatake
++++++++

The datatake ID of source scenes in hexadecimal representation, e.g. 04EBF7.

work_dir
++++++++

``work_dir`` is the main directory in which any subdirectories and files are stored that are generated during processing.
Needs to be provided as full path to an existing directory.

tmp_dir, sar_dir, ard_dir, wbm_dir
++++++++++++++++++++++++++++++++++

Processing creates many intermediate files that are expected to be stored in separate subdirectories. The
default values provided in the example configuration file linked above are recommended and will automatically create
subdirectories relative to the directory specified with ``work_dir``. E.g., ``ard_dir = ARD`` will create the subdirectory
``/<work_dir>/ARD``. Optionally, full paths to existing directories can be provided for all of these parameters.

logfile
+++++++

The path to a log file. If set to ``None``, all logs will be printed to the console.
The file path can be relative to ``work_dir`` or absolute.
Default if not defined: ``None``.

search option I: scene_dir & db_file
++++++++++++++++++++++++++++++++++++

Metadata of any Sentinel-1 scene found in ``scene_dir`` will be stored in an SQLite database file created by :class:`pyrosar.drivers.Archive`.
With ``db_file`` either a full path to an existing database can be provided or it will be created in ``work_dir`` if only
a filename is provided. E.g., ``db_file = scenes.db`` will automatically create the database file ``/<work_dir>/scenes.db``.
``scene_dir`` needs to be provided as full path to an existing directory and will be searched recursively for any Sentinel-1
scenes using the regular expression ``'^S1[AB].*(SAFE|zip)$'``.

search option II: stac_catalog & stac_collections
+++++++++++++++++++++++++++++++++++++++++++++++++

Alternative to searching scenes in a directory and storing their metadata in an SQLite database, scenes can be queried from a STAC catalog.
For this, a STAC URL and one or many collections can be defined with ``stac_catalog`` and ``stac_collections`` respectively.
The scenes are expected to be locally accessible in unpacked folders with the `.SAFE` extension.

kml_file
++++++++

The Sentinel-2 Military Grid Reference System (MGRS) tiling system establishes the basis of the processing chain and a
local reference file containing the respective tile information for processing ARD products is needed. The official
KML file defined for the Sentinel-2 mission provided by ESA can be retrieved `here <https://sentinel.esa.int/documents/247904/1955685/S2A_OPER_GIP_TILPAR_MPC__20151209T095117_V20150622T000000_21000101T000000_B00.kml>`_.
With the ``kml_file`` parameter either a full path to this reference file can be provided or it is expected to be located
in the directory provided with ``work_dir`` if only a filename is provided. E.g., the processor expects to find
``/<work_dir>/s2_grid.kml`` if ``kml_file = s2_grid.kml``.

dem_type
++++++++

Options: ``Copernicus 10m EEA DEM | Copernicus 30m Global DEM II | Copernicus 30m Global DEM | GETASSE30``

The Digital Elevation Model (DEM) that should be used for processing.

Note that water body masks are not available for "GETASSE30", and will therefore not be
included in the product data mask. "Copernicus 10m EEA DEM" and "Copernicus 30m Global DEM II" (both include water body masks)
are retrieved from the `Copernicus Space Component Data Access system (CSCDA) <https://spacedata.copernicus.eu/web/cscda/data-access/registration>`_,
which requires authentication. The processor reads username and password from the environment variables `DEM_USER`
and `DEM_PASS` if possible and otherwise interactively asks for authentication if one of these DEM options is selected.

gdal_threads
++++++++++++

Temporarily changes GDAL_NUM_THREADS during processing. Will be reset after processing has finished.

measurement
+++++++++++

Options: ``gamma | sigma``

The backscatter measurement convention. Either creates gamma naught RTC (:math:`\gamma^0_T`) or sigma naught RTC (:math:`\sigma^0_T`) backscatter.

annotation
++++++++++

A comma-separated list to define the annotation layers to be created for each ARD product.
Supported options:

 + dm: data mask (six masks: not layover not shadow, layover, shadow, ocean, lakes, rivers)
 + ei: ellipsoidal incident angle (needed for computing geolocation accuracy)
 + em: digital elevation model
 + id: acquisition ID image (source scene ID per pixel)
 + lc: RTC local contributing area
 + ld: range look direction angle
 + li: local incident angle
 + np: noise power (NESZ, per polarization)
 + ratio: will automatically be replaced with the following, depending on selected ``measurement``:

   + gs: gamma-sigma ratio: sigma0 RTC / gamma0 RTC (if ``measurement = gamma``)
   + sg: sigma-gamma ratio: gamma0 RTC / sigma0 RTC (if ``measurement = sigma``)

 + wm: wind-modelled backscatter extracted from a Sentinel-1 OCN (ocean) product.
   The sub-product `owiNrcsCmod` is extracted, which is Ocean Wind (OWI) Normalised
   Radar Cross Section (NRCS) predicted using a CMOD model and ECMWF wind model data.
   For each OCN product, a Level-1 counterpart (SLC/GRD) exists.
   The OCN products and corresponding Level-1 products must be searchable in the same way
   via the two search options described above.
   If a sigma naught output layer exists (via ``measurement = sigma`` or `annotation` layer `ratio`),
   a co-polarization wind normalization ratio VRT is created by dividing the measurement by the
   wind-modelled backscatter.

Use one of the following to create no annotation layer:

 + ``annotation =``
 + ``annotation = None``

etad & etad_dir
+++++++++++++++

Determines if the `Extended Timing Annotation Dataset (ETAD) correction <https://sentinel.esa.int/web/sentinel/missions/sentinel-1/data-products/etad-dataset>`_
should be performed or not. If ``etad=True``, ``etad_dir`` is searched for ETAD products matching the respective input SLC
and a new SLC is created in ``tmp_dir``, which is then used for all other processing steps. If ``etad=False``, ``etad_dir``
will be ignored.

Metadata Section
^^^^^^^^^^^^^^^^

format
++++++

A comma-separated list to define the metadata file formats to be created for each ARD product.
Supported options:

 + OGC: XML file according to `OGC EO <https://docs.ogc.org/is/10-157r4/10-157r4.html>`_ standard
 + STAC: JSON file according to the `SpatioTemporal Asset Catalog <https://github.com/radiantearth/stac-spec/>`_ family of specifications

copy_original
+++++++++++++

Copy the original metadata of the source scene(s) into the ARD product directory?
This will copy the manifest.safe file and annotation folder into the subdirectory: ``/source/<ProductIdentifier>``.

access_url, licence, doi & processing_center
++++++++++++++++++++++++++++++++++++++++++++

The metadata files created for each ARD product contain some fields that should not be hidden away and hardcoded with
arbitrary values. Instead, they can be accessed here in order to more easily generate a complete set of metadata. These
fields are mostly relevant if you want to produce ARD products systematically and make them available for others.
If you don't see a need for them you can just leave the fields empty, use the default 'None' or delete this entire section.

Command Line Interface
----------------------
Once a configuration file has been created and all of its parameters have been properly defined, it can be used to start
the processor using the command line interface (CLI) tool provided with the s1ard package.
The tool `s1rb` can be used to create the two radar backscatter products NRB and ORB.

The following options are currently available.

Print a help message for the CLI tool:

::

    s1rb --help

Print the processor version:

::

    s1rb --version

Start the processor using parameters defined in the default section of a ``config.ini`` file:

::

    s1rb -c /path/to/config.ini

Start the processor using parameters defined in section ``SECTION_NAME`` of a ``config.ini`` file:

::

    s1rb -c /path/to/config.ini -s SECTION_NAME

Start the processor using parameters defined in the default section of a ``config.ini`` file but
override some parameters, e.g. ``acq_mode`` and ``annotation``:

::

    s1ard -c /path/to/config.ini --acq_mode IW --annotation dm,id

The argument `snap_gpt_args` is known to require an additional modification so that the `-` characters in the value are not mistaken for argument keys. 
In the example SNAP is instructed to use a maximum of 32GB memory, 20GB cache size and 16 threads.

::

    s1ard -c /path/to/config.ini -- --snap_gpt_args "-J-Xmx32G -c 20G -x -q 16"
