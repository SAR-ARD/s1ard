Folder Structure
================

The following demonstrates a possible structure created to store intermediate and final files during a processor run.
The listed files describe the output of the user configuration parameter ``measurement`` set to ``gamma``
and the following output annotation layers enabled ``annotation = dm, ei, em, id, lc, li, np, ratio``, thus creating S1-NRB products.
The structure is based on the default configuration defined in the `config.ini` file and can be modified by a user.
Folders are highlighted in bold.

.. only:: latex

    This section is currently not supported with LaTeX/PDF as it was written with collapsible elements in HTML.

.. only:: html

    .. collapse:: <b>work_dir</b>
        :open:

            .. collapse:: <b>ard_dir</b>

                .. pull-quote::

                    .. note::

                        The final NRB/ORB tiles sorted into subfolders by MGRS tile.
                        Additional STAC files can be generated using function :func:`s1ard.metadata.stac.make_catalog`:

                        - `collection.json`: a STAC collection file referencing the product-specific STAC item files per MGRS tile
                        - `catalog.json`: a STAC catalog referencing all collections

                    .. collapse:: <b>32TPS</b>

                        .. pull-quote::

                            .. collapse:: <b>S1A_IW_NRB__1SDV_20200103T170705_030639_0382D5_32TPS_8090</b>

                                .. pull-quote::

                                    .. collapse:: <b>annotation</b>

                                        .. pull-quote::

                                            | s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-dm.tif
                                            | s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-ei.tif
                                            | s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-em.tif
                                            | s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-gs.tif
                                            | s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-id.tif
                                            | s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-lc.tif
                                            | s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-li.tif
                                            | s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-np-vh.tif
                                            | s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-np-vv.tif

                                    .. collapse:: <b>measurement</b>

                                        .. pull-quote::

                                            | s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-cc-g-lin.vrt
                                            | s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-vh-g-lin.tif
                                            | s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-vh-g-log.vrt
                                            | s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-vh-s-lin.vrt
                                            | s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-vh-s-log.vrt
                                            | s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-vv-g-lin.tif
                                            | s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-vv-g-log.vrt
                                            | s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-vv-s-lin.vrt
                                            | s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-vv-s-log.vrt

                                    .. collapse:: <b>source</b>

                                        .. pull-quote::

                                            .. note::

                                                This folder contains the metadata of the source product(s). Two files per source
                                                scene contain a subset of the source product metadta in STAC and OGC EO XML format.
                                                The subfolder(s) named after the UID of the source scene(s) further contains all
                                                unaltered metadata found in the source product.

                                            .. collapse:: <b>6A12</b>

                                                .. pull-quote::

                                                    .. collapse:: <b>annotation</b>

                                                        .. pull-quote::

                                                            .. collapse:: <b>calibration</b>

                                                                .. pull-quote::

                                                                    | calibration-s1a-iw1-slc-vh-20200103t170701-20200103t170726-030639-0382d5-001.xml
                                                                    | calibration-s1a-iw1-slc-vv-20200103t170701-20200103t170726-030639-0382d5-004.xml
                                                                    | calibration-s1a-iw2-slc-vh-20200103t170702-20200103t170727-030639-0382d5-002.xml
                                                                    | calibration-s1a-iw2-slc-vv-20200103t170702-20200103t170727-030639-0382d5-005.xml
                                                                    | calibration-s1a-iw3-slc-vh-20200103t170700-20200103t170725-030639-0382d5-003.xml
                                                                    | calibration-s1a-iw3-slc-vv-20200103t170700-20200103t170725-030639-0382d5-006.xml
                                                                    | noise-s1a-iw1-slc-vh-20200103t170701-20200103t170726-030639-0382d5-001.xml
                                                                    | noise-s1a-iw1-slc-vv-20200103t170701-20200103t170726-030639-0382d5-004.xml
                                                                    | noise-s1a-iw2-slc-vh-20200103t170702-20200103t170727-030639-0382d5-002.xml
                                                                    | noise-s1a-iw2-slc-vv-20200103t170702-20200103t170727-030639-0382d5-005.xml
                                                                    | noise-s1a-iw3-slc-vh-20200103t170700-20200103t170725-030639-0382d5-003.xml
                                                                    | noise-s1a-iw3-slc-vv-20200103t170700-20200103t170725-030639-0382d5-006.xml

                                                            | s1a-iw1-slc-vh-20200103t170701-20200103t170726-030639-0382d5-001.xml
                                                            | s1a-iw1-slc-vv-20200103t170701-20200103t170726-030639-0382d5-004.xml
                                                            | s1a-iw2-slc-vh-20200103t170702-20200103t170727-030639-0382d5-002.xml
                                                            | s1a-iw2-slc-vv-20200103t170702-20200103t170727-030639-0382d5-005.xml
                                                            | s1a-iw3-slc-vh-20200103t170700-20200103t170725-030639-0382d5-003.xml
                                                            | s1a-iw3-slc-vv-20200103t170700-20200103t170725-030639-0382d5-006.xml

                                                    | mainfest.safe

                                            | S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12.json
                                            | S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12.xml

                                    .. collapse:: <b>support</b>

                                        .. pull-quote::

                                            | product.xsd
                                            | source.xsd

                                    | S1A_IW_NRB__1SDV_20200103T170705_030639_0382D5_32TPS_8090.json
                                    | S1A_IW_NRB__1SDV_20200103T170705_030639_0382D5_32TPS_8090.xml

                            | ...
                            | collection.json

                    | ...
                    | catalog.json

            .. collapse:: <b>sar_dir</b>

                .. pull-quote::

                    .. note::

                        The SAR processing output and SNAP workflows per source scene.
                        Geocoded products carry an EPSG code suffix.

                    .. collapse:: <b>S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12</b>

                        .. pull-quote::

                            | **S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_geo_32632.data**
                            | S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_geo_32632.dim
                            | S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_geo_32632.xml
                            | S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_gsr.xml
                            | S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_mli.xml
                            | S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_pre.xml
                            | S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_rtc.xml

                            ...

                    ...

            .. collapse:: <b>tmp_dir</b>

                .. pull-quote::

                    .. note::

                        Intermediate non-geocoded SAR processor files per scene.

                        - scene-specific DEM mosaic and intermediate (SNAP) processor files
                        - unpacked ETAD files (\*_ETA_\*)
                        - SLC_etad subfolder: ETAD-corrected SLC

                    .. collapse:: <b>S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12</b>

                        .. pull-quote::

                            | **S1A_IW_ETA__AXDV_20200103T170700_20200103T170727_030639_0382D5_256B.SAFE**
                            | **S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_gsr.data**
                            | **S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_mli.data**
                            | **S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_pre.data**
                            | **S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_rtc.data**
                            | S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_DEM_EEA10.tif
                            | S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_DEM_EEA10.vrt
                            | S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_gsr.dim
                            | S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_gsr.xml
                            | S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_mli.dim
                            | S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_mli.xml
                            | S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_pre.dim
                            | S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_pre.xml
                            | S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_rtc.dim
                            | S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_rtc.xml
                            | ...

                            .. collapse:: <b>SLC_etad</b>

                                .. pull-quote::

                                    **S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12.SAFE**

                    ...

            .. collapse:: <b>wbm_dir</b>

                .. pull-quote::

                    .. note::

                        Water Body Mask tiles in MGRS grid per DEM type.
                        The type/folder names are taken from :func:`pyroSAR.auxdata.dem_autoload`.

                    .. collapse:: <b>Copernicus 10m EEA DEM</b>

                        .. pull-quote::

                            | 32TPR_WBM.tif
                            | 32TPS_WBM.tif
                            | 33TUL_WBM.tif
                            | ...


                    .. collapse:: <b>Copernicus 30m Global DEM II</b>

                        .. pull-quote::

                            | 32TPR_WBM.tif
                            | 32TPS_WBM.tif
                            | 33TUL_WBM.tif
                            | ...

            db_file (optional, if set)
