Changelog
=========

1.5.0 | 2023-10-12
------------------

* Replace `gs` and `sg` annotation options with `ratio` (`#116 <https://github.com/SAR-ARD/S1_NRB/pull/116>`_)
* Metadata/review (`#117 <https://github.com/SAR-ARD/S1_NRB/pull/117>`_)
* Equivalent Number of Looks (`#113 <https://github.com/SAR-ARD/S1_NRB/pull/113>`_)
* [copy_src_meta] fixed bug in reading zip content on Windows (`#124 <https://github.com/SAR-ARD/S1_NRB/pull/124>`_)
* Documentation: Table of abbreviations (`#123 <https://github.com/SAR-ARD/S1_NRB/pull/123>`_)
* fixed bug in GRD buffering of ascending scenes (`#126 <https://github.com/SAR-ARD/S1_NRB/pull/126>`_)
* new annotation layer "range look direction angle" (`#103 <https://github.com/SAR-ARD/S1_NRB/pull/103>`_)
* ENL calculation: Suppress warnings and increase default block_size (`#127 <https://github.com/SAR-ARD/S1_NRB/pull/127>`_)
* Add missing pyproj dependency (`#128 <https://github.com/SAR-ARD/S1_NRB/pull/128>`_)
* Simplified datamask for ORB product (`#122 <https://github.com/SAR-ARD/S1_NRB/pull/122>`_)
* Update .readthedocs.yaml (`#129 <https://github.com/SAR-ARD/S1_NRB/pull/129>`_)
* [nrb.create_vrt] fixed bug in handling default 'options=None' (`#132 <https://github.com/SAR-ARD/S1_NRB/pull/132>`_)
* [docs] point to right environment.yaml when installing specific version (`#133 <https://github.com/SAR-ARD/S1_NRB/pull/133>`_)
* Fix missing STAC FileExtension entries (`#131 <https://github.com/SAR-ARD/S1_NRB/pull/131>`_)
* Accommodate ORB product (`#121 <https://github.com/SAR-ARD/S1_NRB/pull/121>`_)
* rename config default annotation IDs gs and sg to ratio (`#135 <https://github.com/SAR-ARD/S1_NRB/pull/135>`_)
* [snap.process] skip GRD buffering if list is empty (`#139 <https://github.com/SAR-ARD/S1_NRB/pull/139>`_)
* Refer to original source metadata in source XML and JSON (`#136 <https://github.com/SAR-ARD/S1_NRB/pull/136>`_)
* wind normalization (`#138 <https://github.com/SAR-ARD/S1_NRB/pull/138>`_)
* Look direction angle improvements (`#141 <https://github.com/SAR-ARD/S1_NRB/pull/141>`_)
* do not look for source metadata files if copying is not user-configured (`#142 <https://github.com/SAR-ARD/S1_NRB/pull/142>`_)
* change EW spacing from 20 to 40 m (`#143 <https://github.com/SAR-ARD/S1_NRB/pull/143>`_)
* XML product metadata improvements (`#137 <https://github.com/SAR-ARD/S1_NRB/pull/137>`_)
* Metadata/review (`#140 <https://github.com/SAR-ARD/S1_NRB/pull/140>`_)
* wind normalization - GDAL options (`#144 <https://github.com/SAR-ARD/S1_NRB/pull/144>`_)
* Require pyroSAR >=0.22.0 and update license year (`#145 <https://github.com/SAR-ARD/S1_NRB/pull/145>`_)
* documentation improvements (`#146 <https://github.com/SAR-ARD/S1_NRB/pull/146>`_)
* STACArchive file path handling (`#148 <https://github.com/SAR-ARD/S1_NRB/pull/148>`_)
* geometry buffering for minimum overlap (`#147 <https://github.com/SAR-ARD/S1_NRB/pull/147>`_)
* apply RTC to sigma0 (`#149 <https://github.com/SAR-ARD/S1_NRB/pull/149>`_)
* config 'mode': removed 'all', added 'orb'; renamed module 'nrb' to 'ard' (`#150 <https://github.com/SAR-ARD/S1_NRB/pull/150>`_)

`Full v1.5.0 Changelog <https://github.com/SAR-ARD/S1_NRB/compare/v1.4.0...v1.5.0>`_

1.4.0 | 2023-07-04
------------------

* various bug fixes (`#94 <https://github.com/SAR-ARD/S1_NRB/pull/94>`_)
* datatake gap handling (`#95 <https://github.com/SAR-ARD/S1_NRB/pull/95>`_)
* new configuration parameter 'datatake' (`#96 <https://github.com/SAR-ARD/S1_NRB/pull/96>`_)
* increased STAC access robustness (`#97 <https://github.com/SAR-ARD/S1_NRB/pull/97>`_)
* STACArchive bug fixes (`#98 <https://github.com/SAR-ARD/S1_NRB/pull/98>`_)
* Optional `datatake` parameter (`#99 <https://github.com/SAR-ARD/S1_NRB/pull/99>`_)
* bug fixes (`#100 <https://github.com/SAR-ARD/S1_NRB/pull/100>`_)
* Bug fix to allow `annotation = None` (`#102 <https://github.com/SAR-ARD/S1_NRB/pull/102>`_)
* Save original source metadata  (`#104 <https://github.com/SAR-ARD/S1_NRB/pull/104>`_)
* do not continue on error (`#105 <https://github.com/SAR-ARD/S1_NRB/pull/105>`_)
* Always use ESA border noise removal (`#106 <https://github.com/SAR-ARD/S1_NRB/pull/106>`_)
* [nrb] remove dataset if mask is nodata-only (`#108 <https://github.com/SAR-ARD/S1_NRB/pull/108>`_)
* Bug fix: Save original source metadata (`#109 <https://github.com/SAR-ARD/S1_NRB/pull/109>`_)
* New metadata config parameters (`#110 <https://github.com/SAR-ARD/S1_NRB/pull/110>`_)
* support for scenes acquired in NRT Slicing mode (`#112 <https://github.com/SAR-ARD/S1_NRB/pull/112>`_)

`Full Changelog <https://github.com/SAR-ARD/S1_NRB/compare/v1.3.0...v1.4.0>`_

1.3.0 | 2023-05-24
------------------

* SNAP RTC: increase DEM oversampling by a factor of two (`#78 <https://github.com/SAR-ARD/S1_NRB/pull/78>`_)
* nrb.format: do not hardcode src_nodata and read it from the data instead (`#79 <https://github.com/SAR-ARD/S1_NRB/pull/79>`_)
* enable configuration via command line arguments (`#80 <https://github.com/SAR-ARD/S1_NRB/pull/80>`_)
* improved date parsing (`#81 <https://github.com/SAR-ARD/S1_NRB/pull/81>`_)
* scene search via STAC (`#82 <https://github.com/SAR-ARD/S1_NRB/pull/82>`_)
* enhanced time filtering (`#84 <https://github.com/SAR-ARD/S1_NRB/pull/84>`_)
* general processor improvements (`#85 <https://github.com/SAR-ARD/S1_NRB/pull/85>`_)

`Full Changelog <https://github.com/SAR-ARD/S1_NRB/compare/v1.2.0...v1.3.0>`_

1.2.0 | 2022-12-29
------------------

* improved geometry handling (`#71 <https://github.com/SAR-ARD/S1_NRB/pull/71>`_)
* DEM handling improvements (`#72 <https://github.com/SAR-ARD/S1_NRB/pull/72>`_)
* GRD buffering by (`#73 <https://github.com/SAR-ARD/S1_NRB/pull/73>`_)
* add DEM as additional output layer (`#70 <https://github.com/SAR-ARD/S1_NRB/pull/70>`_)
* sigma0 processing and annotation layer configuration (`#74 <https://github.com/SAR-ARD/S1_NRB/pull/74>`_)

`Full Changelog <https://github.com/SAR-ARD/S1_NRB/compare/v1.1.0...v1.2.0>`_

1.1.0 | 2022-09-29
------------------

* documentation improvements (`#60 <https://github.com/SAR-ARD/S1_NRB/pull/60>`_)
* installation update (`#61 <https://github.com/SAR-ARD/S1_NRB/pull/61>`_)
* Process restructuring (`#63 <https://github.com/SAR-ARD/S1_NRB/pull/63>`_)
* minor structural changes and bug fixes (`#65 <https://github.com/SAR-ARD/S1_NRB/pull/65>`_)
* documentation update reflecting the recent process restructuring (`#66 <https://github.com/SAR-ARD/S1_NRB/pull/66>`_)
* renamed processing mode 'snap' to 'rtc' (`#67 <https://github.com/SAR-ARD/S1_NRB/pull/67>`_)

`Full Changelog <https://github.com/SAR-ARD/S1_NRB/compare/v1.0.2...v1.1.0>`_

1.0.2 | 2022-08-24
------------------

* Fix error in handling of temporary VRTs (`#50 <https://github.com/SAR-ARD/S1_NRB/pull/50>`_)
* Adjustments to VRT log scaling (`#52 <https://github.com/SAR-ARD/S1_NRB/pull/52>`_)
* [metadata] read nodata values directly from files (instead of hard-coding them) (`#53 <https://github.com/SAR-ARD/S1_NRB/pull/53>`_)
* use type identifier in scene-specific DEM file names (`#55 <https://github.com/SAR-ARD/S1_NRB/pull/55>`_)
* Add VRT assets to STAC files (`#56 <https://github.com/SAR-ARD/S1_NRB/pull/56>`_)
* Fix and improve metadata geometry handling (`#57 <https://github.com/SAR-ARD/S1_NRB/pull/57>`_)
* SNAP 9 compatibility (`#58 <https://github.com/SAR-ARD/S1_NRB/pull/58>`_)

`Full Changelog <https://github.com/SAR-ARD/S1_NRB/compare/v1.0.1...v1.0.2>`_

1.0.1 | 2022-07-03
------------------

* dem handling improvements (`#45 <https://github.com/SAR-ARD/S1_NRB/pull/45>`_)

`Full Changelog <https://github.com/SAR-ARD/S1_NRB/compare/v1.0.0...v1.0.1>`_

1.0.0 | 2022-06-23
------------------

* Dockerfile to build S1_NRB image (`#31 <https://github.com/SAR-ARD/S1_NRB/pull/31>`_)
* adjustments to nodata value (`#28 <https://github.com/SAR-ARD/S1_NRB/pull/28>`_)
* renamed XML tag 'nrb' to 's1-nrb' (`#36 <https://github.com/SAR-ARD/S1_NRB/pull/36>`_)
* Metadata & Config Improvements (`#30 <https://github.com/SAR-ARD/S1_NRB/pull/30>`_)
* Geolocation accuracy (`#40 <https://github.com/SAR-ARD/S1_NRB/pull/40>`_)
* various bug fixes and documentation improvements

`Full Changelog <https://github.com/SAR-ARD/S1_NRB/compare/v0.4.2...v1.0.0>`_

0.4.2 | 2022-06-16
------------------

* Update documentation (`#27 <https://github.com/SAR-ARD/S1_NRB/pull/27>`_)
* find unpacked .SAFE scenes in scene_dir (instead of just .zip) (`aea53a5 <https://github.com/SAR-ARD/S1_NRB/commit/aea53a57bc5fa1418fea4f46f69b41b7332909b1>`_)

`Full Changelog <https://github.com/SAR-ARD/S1_NRB/compare/v0.4.1...v0.4.2>`_

0.4.1 | 2022-06-01
------------------

* handle ETAD products as zip, tar, and SAFE (`#25 <https://github.com/SAR-ARD/S1_NRB/pull/25>`_)
* set dem download authentication via env. variables (`#26 <https://github.com/SAR-ARD/S1_NRB/pull/26>`_)
* various bug fixes

`Full Changelog <https://github.com/SAR-ARD/S1_NRB/compare/v0.4.0...v0.4.1>`_

0.4.0 | 2022-05-30
------------------

* outsourced and restructured DEM preparation functionality (`#18 <https://github.com/SAR-ARD/S1_NRB/pull/18>`_)
* outsourced ETAD correction to dedicated module (`#19 <https://github.com/SAR-ARD/S1_NRB/pull/19>`_)
* XML validation & improvements (`#17 <https://github.com/SAR-ARD/S1_NRB/pull/17>`_)
* Restructuring and cleanup (`#20 <https://github.com/SAR-ARD/S1_NRB/pull/20>`_)
* outsourced NRB formatting to dedicated module (`#21 <https://github.com/SAR-ARD/S1_NRB/pull/21>`_)
* extended acquisition mode support (`#22 <https://github.com/SAR-ARD/S1_NRB/pull/22>`_)
* Set up sphinx documentation (`#23 <https://github.com/SAR-ARD/S1_NRB/pull/23>`_)
* AOI scene selection (`#24 <https://github.com/SAR-ARD/S1_NRB/pull/24>`_)

`Full Changelog <https://github.com/SAR-ARD/S1_NRB/compare/v0.3.0...v0.4.0>`_

0.3.0 | 2022-03-30
------------------

* Updated metadata module (`#9 <https://github.com/SAR-ARD/S1_NRB/pull/9>`_)
* Modified `prepare_dem` interface (`#10 <https://github.com/SAR-ARD/S1_NRB/pull/10>`_)
* Various improvements (`#11 <https://github.com/SAR-ARD/S1_NRB/pull/11>`_)
* Modified working directory structure (`#12 <https://github.com/SAR-ARD/S1_NRB/pull/12>`_)
* Updated `ancillary.py` (`#13 <https://github.com/SAR-ARD/S1_NRB/pull/13>`_)
* Added ETAD correction (`#14 <https://github.com/SAR-ARD/S1_NRB/pull/14>`_)
* Improved RGB composite (`#15 <https://github.com/SAR-ARD/S1_NRB/pull/15>`_)
* Store DEM/WBM tiles in UTM zones different to the native MGRS zone (`#16 <https://github.com/SAR-ARD/S1_NRB/pull/16>`_)

`Full Changelog <https://github.com/SAR-ARD/S1_NRB/compare/v0.2.0...v0.3.0>`_

0.2.0 | 2022-03-03
------------------

`Full Changelog <https://github.com/SAR-ARD/S1_NRB/compare/v0.1.0...v0.2.0>`_

0.1.0 | 2022-01-14
------------------

