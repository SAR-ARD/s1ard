Changelog
=========

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

