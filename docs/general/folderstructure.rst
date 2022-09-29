Folder Structure
================

The following describes the structure created to store intermediate and final files during a processor run.
The structure is based on the default configuration defined in the `config.ini` file and can be modified by a user.
Folders are highlighted in bold.

.. only:: latex

    This section is currently not supported with LaTeX/PDF as it was written with collapsible elements in HTML.

.. only:: html

    .. collapse:: <b>work_dir</b>
        :open:

            .. collapse:: <b>log_dir</b>

                .. pull-quote::

                    .. note::

                        Log files for each processor run containing the full processor configuration (`config.ini`),
                        the versions of relevant installed software, and details on individual processing steps.

                    | 20220719T1339_process.log
                    | 20220719T1032_process.log
                    | 20220708T1118_process.log
                    | ...

            .. collapse:: <b>nrb_dir</b>

                .. pull-quote::

                    .. note::

                        The final S1-NRB tiles sorted into subfolders by MGRS tile.
                        Additional STAC files can be generated using function :func:`S1_NRB.metadata.stac.make_catalog`:

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
                                            | s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-vh-g-log.tif
                                            | s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-vh-s-lin.tif
                                            | s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-vh-s-log.tif
                                            | s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-vv-g-lin.tif
                                            | s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-vv-g-log.tif
                                            | s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-vv-s-lin.tif
                                            | s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-vv-s-log.tif

                                    .. collapse:: <b>source</b>

                                        .. pull-quote::

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

            .. collapse:: <b>rtc_dir</b>

                .. pull-quote::

                    .. note::

                        The RTC processing output and SNAP workflows per source scene.
                        Geocoded products carry an EPSG code suffix.

                    .. collapse:: <b>S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12</b>

                        .. pull-quote::

                            | S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_geo_32632.data
                            | S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_geo_32632.dim
                            | S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_geo_32632.xml
                            | S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_gsr.xml
                            | S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_mli.xml
                            | S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_rtc.xml

                            ...

                    ...

            .. collapse:: <b>tmp_dir</b>

                .. pull-quote::

                    .. note::

                        Intermediate non-geocoded RTC processor files per scene.

                        - scene-specific DEM mosaic and intermediate (SNAP) processor files
                        - unpacked ETAD files (\*_ETA_\*)
                        - SLC_etad subfolder: ETAD-corrected SLCs

                    .. collapse:: <b>S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12</b>

                        .. pull-quote::

                            | S1A_IW_ETA__AXDV_20200103T170700_20200103T170727_030639_0382D5_256B.SAFE
                            | S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_gsr.data
                            | S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_mli.data
                            | S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_rtc.data
                            | S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_DEM_EEA10.tif
                            | S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_DEM_EEA10.vrt
                            | S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_gsr.dim
                            | S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_gsr.xml
                            | S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_mli.dim
                            | S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_mli.xml
                            | S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_rtc.dim
                            | S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_rtc.xml
                            | ...

                            .. collapse:: <b>SLC_etad</b>

                                .. pull-quote::

                                    S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12.SAFE

                    ...

            .. collapse:: <b>wbm_dir</b>

                .. pull-quote::

                    .. note::

                        Water Body Mask tiles in MGRS grid per DEM type.
                        The type names are taken from :func:`pyroSAR.auxdata.dem_autoload`.

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

            scenes.db
