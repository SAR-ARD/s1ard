Changelog
=========

2.6.2 | 2025-12-08
------------------

* typing: do not compare lxml.ElementTree using | but with Union (`#281 <https://github.com/SAR-ARD/s1ard/pull/281>`_)


`Full v2.6.2 Changelog <https://github.com/SAR-ARD/s1ard/compare/v2.6.1...v2.6.2>`_

2.6.1 | 2025-11-25
------------------

* ard.product_info: bug fix, renamed par dir_out to dir_ard (`#279 <https://github.com/SAR-ARD/s1ard/pull/279>`_)


`Full v2.6.1 Changelog <https://github.com/SAR-ARD/s1ard/compare/v2.6.0...v2.6.1>`_

2.6.0 | 2025-11-24
------------------

* added python-publish GitHub Actions workflow (`#275 <https://github.com/SAR-ARD/s1ard/pull/275>`_)
* import from cesard (`#276 <https://github.com/SAR-ARD/s1ard/pull/276>`_)


`Full v2.6.0 Changelog <https://github.com/SAR-ARD/s1ard/compare/v2.5.0...v2.6.0>`_

2.5.0 | 2025-10-27
------------------

* `dem.prepare`: reduce UTM DEM extent (`#273 <https://github.com/SAR-ARD/s1ard/pull/273>`_)
* metadata generalization (`#272 <https://github.com/SAR-ARD/s1ard/pull/272>`_)


`Full v2.5.0 Changelog <https://github.com/SAR-ARD/s1ard/compare/v2.4.0...v2.5.0>`_

2.4.0 | 2025-09-23
------------------

* dem interface modifications (`#263 <https://github.com/SAR-ARD/s1ard/pull/263>`_)
* MGRS grid KML: changed URL and download procedure (`#270 <https://github.com/SAR-ARD/s1ard/pull/270>`_)
* processor modularization, DEM improvements, spacing checks (`#261 <https://github.com/SAR-ARD/s1ard/pull/261>`_)

`Full v2.4.0 Changelog <https://github.com/SAR-ARD/s1ard/compare/v2.3.2...v2.4.0>`_

2.3.2 | 2025-08-01
------------------

* `calc_product_start_stop`: fixed bug in values exceeding valid range (`#267 <https://github.com/SAR-ARD/s1ard/pull/267>`_)
* `snap.pre`: remove grd border boise for all sentinel-1 satellites (`#268 <https://github.com/SAR-ARD/s1ard/pull/268>`_)

`Full v2.3.2 Changelog <https://github.com/SAR-ARD/s1ard/compare/v2.3.1...v2.3.2>`_

2.3.1 | 2025-07-07
------------------

* `metadata.extract`: OSV info gathering improvements (`#264 <https://github.com/SAR-ARD/s1ard/pull/264>`_)
* `ard.calc_product_start_stop`: remove duplicate grid points (`#265 <https://github.com/SAR-ARD/s1ard/pull/265>`_)

`Full v2.3.1 Changelog <https://github.com/SAR-ARD/s1ard/compare/v2.3.0...v2.3.1>`_

2.3.0 | 2025-05-26
------------------

* `snap.process`: fixed bug in not reassigning `out_pre` when buffering has already been performed (`#245 <https://github.com/SAR-ARD/s1ard/pull/245>`_)
* STAC: do not insert links for radiometric and geometric accuracy reference if the value is None (`#246 <https://github.com/SAR-ARD/s1ard/pull/246>`_)
* STAC: do not create link fields when URL is None (`#247 <https://github.com/SAR-ARD/s1ard/pull/247>`_)
* require at least one pixel overlap with the MGRS tile (`#248 <https://github.com/SAR-ARD/s1ard/pull/248>`_)
* product acquisition time improvements (`#251 <https://github.com/SAR-ARD/s1ard/pull/251>`_)
* `ard.calc_product_start_stop`: enable extrapolation (`#252 <https://github.com/SAR-ARD/s1ard/pull/252>`_)
* processing per datatake, slice number filtering (`#250 <https://github.com/SAR-ARD/s1ard/pull/250>`_)
* `ard.calc_product_start_stop`: use geopandas for intersection (`#255 <https://github.com/SAR-ARD/s1ard/pull/255>`_)
* scene and tile search improvements (`#254 <https://github.com/SAR-ARD/s1ard/pull/254>`_)
* support scenes in individual composition (`#256 <https://github.com/SAR-ARD/s1ard/pull/256>`_)
* support for Sentinel-1C (and -1D) (`#257 <https://github.com/SAR-ARD/s1ard/pull/257>`_)
* `config`: new function `init`, write bug fix (`#259 <https://github.com/SAR-ARD/s1ard/pull/259>`_)
* `search.check_acquisition_completeness`: bug fix (`#260 <https://github.com/SAR-ARD/s1ard/pull/260>`_)
* metadata extraction: outsourced URLs to mapping, changed to https for some (`#249 <https://github.com/SAR-ARD/s1ard/pull/249>`_)

`Full v2.3.0 Changelog <https://github.com/SAR-ARD/s1ard/compare/v2.2.0...v2.3.0>`_

2.2.0 | 2025-02-21
------------------

* scene search optimizations (`#241 <https://github.com/SAR-ARD/s1ard/pull/241>`_)

`Full v2.2.0 Changelog <https://github.com/SAR-ARD/s1ard/compare/v2.1.0...v2.2.0>`_

2.1.0 | 2025-01-21
------------------

* documentation: usage description improvements (`#209 <https://github.com/SAR-ARD/s1ard/pull/209>`_)
* add file hash to STAC metadata (`#210 <https://github.com/SAR-ARD/s1ard/pull/210>`_)
* remove config `kml_file` and auto-download with `ancillary.get_kml` (`#211 <https://github.com/SAR-ARD/s1ard/pull/211>`_)
* `ancillary.get_max_ext`: check geometries' CRSs are the same, new arg `crs` (`#212 <https://github.com/SAR-ARD/s1ard/pull/212>`_)
* make config parameter `scene_dir` optional (`#213 <https://github.com/SAR-ARD/s1ard/pull/213>`_)
* fix readthedocs build (`#216 <https://github.com/SAR-ARD/s1ard/pull/216>`_)
* improve geometry handling in search queries (`#214 <https://github.com/SAR-ARD/s1ard/pull/214>`_)
* start processing with `s1rb process`; new option `s1rb init` to initialize config file (`#215 <https://github.com/SAR-ARD/s1ard/pull/215>`_)
* require `numpy<2.0` (`#219 <https://github.com/SAR-ARD/s1ard/pull/219>`_)
* introduce geoparquet scene search (`#217 <https://github.com/SAR-ARD/s1ard/pull/217>`_)
* `metadata.extract.calc_enl`: raise error when array is empty (`#221 <https://github.com/SAR-ARD/s1ard/pull/221>`_)
* use COG media type for VRTs (`#222 <https://github.com/SAR-ARD/s1ard/pull/222>`_)
* `snap`: skip GRD buffering if slice number cannot be determined (`#225 <https://github.com/SAR-ARD/s1ard/pull/225>`_)
* limit pystac to <1.11.0 (`#223 <https://github.com/SAR-ARD/s1ard/pull/223>`_)
* valid data mask improvements (`#226 <https://github.com/SAR-ARD/s1ard/pull/226>`_)
* improved file locking during SNAP processing (`#227 <https://github.com/SAR-ARD/s1ard/pull/227>`_)
* `metadata.extract.calc_enl`: return None if ENL cannot be computed (`#228 <https://github.com/SAR-ARD/s1ard/pull/228>`_)
* improved tmp file handling (`#230 <https://github.com/SAR-ARD/s1ard/pull/230>`_)
* `check_acquisition_completeness`: catch case when ASF result is empty (`#229 <https://github.com/SAR-ARD/s1ard/pull/229>`_)
* docs: corrected installation path, added missing links (`#232 <https://github.com/SAR-ARD/s1ard/pull/232>`_)
* file locking improvements (`#234 <https://github.com/SAR-ARD/s1ard/pull/234>`_)
* new function `tile_extraction.wkt_to_geom` (`#237 <https://github.com/SAR-ARD/s1ard/pull/237>`_)
* `extract._vec_from_srccoords`: code improvements (`#238 <https://github.com/SAR-ARD/s1ard/pull/238>`_)
* `search.collect_neighbors`: filter if more than two neighbors are found (`#236 <https://github.com/SAR-ARD/s1ard/pull/236>`_)
* snap pre modifications (`#239 <https://github.com/SAR-ARD/s1ard/pull/239>`_)

`Full v2.1.0 Changelog <https://github.com/SAR-ARD/s1ard/compare/v2.0.0...v2.1.0>`_

2.0.0 | 2024-05-24
------------------

* avoid permission error when passing tempfile to gdalwarp (`#171 <https://github.com/SAR-ARD/s1ard/pull/171>`_)
* layover-shadow mask: handle value 3 (layover and shadow) (`#177 <https://github.com/SAR-ARD/s1ard/pull/177>`_)
* new dedicated search module (`#172 <https://github.com/SAR-ARD/s1ard/pull/172>`_)
* ard.create_vrt: new argument 'dtype' (`#179 <https://github.com/SAR-ARD/s1ard/pull/179>`_)
* [search.scene_select] parse date only if necessary (`#180 <https://github.com/SAR-ARD/s1ard/pull/180>`_)
* [search.scene_select] compare strings with '==' instead of 'is' (`#181 <https://github.com/SAR-ARD/s1ard/pull/181>`_)
* [search.scene_select] avoid multiple geometry definitions (`#183 <https://github.com/SAR-ARD/s1ard/pull/183>`_)
* fixed handling of empty list of neighbors during SLC processing (`#186 <https://github.com/SAR-ARD/s1ard/pull/186>`_)
* [search.scene_select] fixed incomplete result for 'aoi_geometry' and 'vectorobject' (`#188 <https://github.com/SAR-ARD/s1ard/pull/188>`_)
* [aoi_from_scene] fixed bug in handling scenes crossing the equator (`#187 <https://github.com/SAR-ARD/s1ard/pull/187>`_)
* STAC query optimizations (`#185 <https://github.com/SAR-ARD/s1ard/pull/185>`_)
* [processor] fixed 'empty selection' message (`#190 <https://github.com/SAR-ARD/s1ard/pull/190>`_)
* new configuration parameter 'scene' (`#184 <https://github.com/SAR-ARD/s1ard/pull/184>`_)
* scene search bug fixes (`#191 <https://github.com/SAR-ARD/s1ard/pull/191>`_)
* use file locking for SNAP processing (`#192 <https://github.com/SAR-ARD/s1ard/pull/192>`_)
* modernize build process (`#194 <https://github.com/SAR-ARD/s1ard/pull/194>`_)
* fix deprecated conda --force in docker (`#193 <https://github.com/SAR-ARD/s1ard/pull/193>`_)
* support for SNAP 10 (`#195 <https://github.com/SAR-ARD/s1ard/pull/195>`_)
* added dedicated documentation section on scene search (`#196 <https://github.com/SAR-ARD/s1ard/pull/196>`_)
* [search.asf_select] ensure naive datetime objects are defined as UTC (`#197 <https://github.com/SAR-ARD/s1ard/pull/197>`_)
* [search.ASF.scanMetadata] date formatting bug fix (`#198 <https://github.com/SAR-ARD/s1ard/pull/198>`_)
* rename package (`#199 <https://github.com/SAR-ARD/s1ard/pull/199>`_)
* update documentation links (`#200 <https://github.com/SAR-ARD/s1ard/pull/200>`_)
* replaced configuration `log_dir` with `logfile`, cleaned up logging (`#201 <https://github.com/SAR-ARD/s1ard/pull/201>`_)
* renamed the command line tool from s1ard to s1rb (`#202 <https://github.com/SAR-ARD/s1ard/pull/202>`_)

`Full v2.0.0 Changelog <https://github.com/SAR-ARD/s1ard/compare/v1.6.2...v2.0.0>`_

1.6.2 | 2023-11-23
------------------

* Update metadata links (`#165 <https://github.com/SAR-ARD/s1ard/pull/165>`_)
* Fix missing datamask layers in metadata (`#164 <https://github.com/SAR-ARD/s1ard/pull/164>`_)
* Add wind normalisation metadata fields (`#166 <https://github.com/SAR-ARD/s1ard/pull/166>`_)
* documentation updates (`#167 <https://github.com/SAR-ARD/s1ard/pull/167>`_)
* [metadata.xml.product_xml] add geo acc. reference only if performed (`#168 <https://github.com/SAR-ARD/s1ard/pull/168>`_)
* require pyroSAR>=0.23.0 (`#169 <https://github.com/SAR-ARD/s1ard/pull/169>`_)


`Full v1.6.2 Changelog <https://github.com/SAR-ARD/s1ard/compare/v1.6.1...v1.6.2>`_

1.6.1 | 2023-11-17
------------------

* use relative paths in wind normalization VRT (`#163 <https://github.com/SAR-ARD/s1ard/pull/163>`_)

`Full v1.6.1 Changelog <https://github.com/SAR-ARD/s1ard/compare/v1.6.0...v1.6.1>`_

1.6.0 | 2023-11-15
------------------

* central documentation literature management (`#151 <https://github.com/SAR-ARD/s1ard/pull/151>`_)
* Use the official Continuum Docker base image (`#152 <https://github.com/SAR-ARD/s1ard/pull/152>`_)
* re-introduce recently lost radiometric terrain correction (`#154 <https://github.com/SAR-ARD/s1ard/pull/154>`_)
* strip line breaks from all parameters passed via the command line (`#155 <https://github.com/SAR-ARD/s1ard/pull/155>`_)
* increase OCN gap fill distance (`#156 <https://github.com/SAR-ARD/s1ard/pull/156>`_)
* data mask modifications (`#157 <https://github.com/SAR-ARD/s1ard/pull/157>`_)
* [config] corrected list of allowed modes (`#158 <https://github.com/SAR-ARD/s1ard/pull/158>`_)
* search OCN scenes by buffered start and stop time (`#160 <https://github.com/SAR-ARD/s1ard/pull/160>`_)
* separate ocean, rivers and lakes into separate data mask bands (`#161 <https://github.com/SAR-ARD/s1ard/pull/161>`_)

`Full v1.6.0 Changelog <https://github.com/SAR-ARD/s1ard/compare/v1.5.0...v1.6.0>`_

1.5.0 | 2023-10-12
------------------

* Replace `gs` and `sg` annotation options with `ratio` (`#116 <https://github.com/SAR-ARD/s1ard/pull/116>`_)
* Metadata/review (`#117 <https://github.com/SAR-ARD/s1ard/pull/117>`_)
* Equivalent Number of Looks (`#113 <https://github.com/SAR-ARD/s1ard/pull/113>`_)
* [copy_src_meta] fixed bug in reading zip content on Windows (`#124 <https://github.com/SAR-ARD/s1ard/pull/124>`_)
* Documentation: Table of abbreviations (`#123 <https://github.com/SAR-ARD/s1ard/pull/123>`_)
* fixed bug in GRD buffering of ascending scenes (`#126 <https://github.com/SAR-ARD/s1ard/pull/126>`_)
* new annotation layer "range look direction angle" (`#103 <https://github.com/SAR-ARD/s1ard/pull/103>`_)
* ENL calculation: Suppress warnings and increase default block_size (`#127 <https://github.com/SAR-ARD/s1ard/pull/127>`_)
* Add missing pyproj dependency (`#128 <https://github.com/SAR-ARD/s1ard/pull/128>`_)
* Simplified datamask for ORB product (`#122 <https://github.com/SAR-ARD/s1ard/pull/122>`_)
* Update .readthedocs.yaml (`#129 <https://github.com/SAR-ARD/s1ard/pull/129>`_)
* [nrb.create_vrt] fixed bug in handling default 'options=None' (`#132 <https://github.com/SAR-ARD/s1ard/pull/132>`_)
* [docs] point to right environment.yaml when installing specific version (`#133 <https://github.com/SAR-ARD/s1ard/pull/133>`_)
* Fix missing STAC FileExtension entries (`#131 <https://github.com/SAR-ARD/s1ard/pull/131>`_)
* Accommodate ORB product (`#121 <https://github.com/SAR-ARD/s1ard/pull/121>`_)
* rename config default annotation IDs gs and sg to ratio (`#135 <https://github.com/SAR-ARD/s1ard/pull/135>`_)
* [snap.process] skip GRD buffering if list is empty (`#139 <https://github.com/SAR-ARD/s1ard/pull/139>`_)
* Refer to original source metadata in source XML and JSON (`#136 <https://github.com/SAR-ARD/s1ard/pull/136>`_)
* wind normalization (`#138 <https://github.com/SAR-ARD/s1ard/pull/138>`_)
* Look direction angle improvements (`#141 <https://github.com/SAR-ARD/s1ard/pull/141>`_)
* do not look for source metadata files if copying is not user-configured (`#142 <https://github.com/SAR-ARD/s1ard/pull/142>`_)
* change EW spacing from 20 to 40 m (`#143 <https://github.com/SAR-ARD/s1ard/pull/143>`_)
* XML product metadata improvements (`#137 <https://github.com/SAR-ARD/s1ard/pull/137>`_)
* Metadata/review (`#140 <https://github.com/SAR-ARD/s1ard/pull/140>`_)
* wind normalization - GDAL options (`#144 <https://github.com/SAR-ARD/s1ard/pull/144>`_)
* Require pyroSAR >=0.22.0 and update license year (`#145 <https://github.com/SAR-ARD/s1ard/pull/145>`_)
* documentation improvements (`#146 <https://github.com/SAR-ARD/s1ard/pull/146>`_)
* STACArchive file path handling (`#148 <https://github.com/SAR-ARD/s1ard/pull/148>`_)
* geometry buffering for minimum overlap (`#147 <https://github.com/SAR-ARD/s1ard/pull/147>`_)
* apply RTC to sigma0 (`#149 <https://github.com/SAR-ARD/s1ard/pull/149>`_)
* config 'mode': removed 'all', added 'orb'; renamed module 'nrb' to 'ard' (`#150 <https://github.com/SAR-ARD/s1ard/pull/150>`_)

`Full v1.5.0 Changelog <https://github.com/SAR-ARD/s1ard/compare/v1.4.0...v1.5.0>`_

1.4.0 | 2023-07-04
------------------

* various bug fixes (`#94 <https://github.com/SAR-ARD/s1ard/pull/94>`_)
* datatake gap handling (`#95 <https://github.com/SAR-ARD/s1ard/pull/95>`_)
* new configuration parameter 'datatake' (`#96 <https://github.com/SAR-ARD/s1ard/pull/96>`_)
* increased STAC access robustness (`#97 <https://github.com/SAR-ARD/s1ard/pull/97>`_)
* STACArchive bug fixes (`#98 <https://github.com/SAR-ARD/s1ard/pull/98>`_)
* Optional `datatake` parameter (`#99 <https://github.com/SAR-ARD/s1ard/pull/99>`_)
* bug fixes (`#100 <https://github.com/SAR-ARD/s1ard/pull/100>`_)
* Bug fix to allow `annotation = None` (`#102 <https://github.com/SAR-ARD/s1ard/pull/102>`_)
* Save original source metadata  (`#104 <https://github.com/SAR-ARD/s1ard/pull/104>`_)
* do not continue on error (`#105 <https://github.com/SAR-ARD/s1ard/pull/105>`_)
* Always use ESA border noise removal (`#106 <https://github.com/SAR-ARD/s1ard/pull/106>`_)
* [nrb] remove dataset if mask is nodata-only (`#108 <https://github.com/SAR-ARD/s1ard/pull/108>`_)
* Bug fix: Save original source metadata (`#109 <https://github.com/SAR-ARD/s1ard/pull/109>`_)
* New metadata config parameters (`#110 <https://github.com/SAR-ARD/s1ard/pull/110>`_)
* support for scenes acquired in NRT Slicing mode (`#112 <https://github.com/SAR-ARD/s1ard/pull/112>`_)

`Full v1.4.0 Changelog <https://github.com/SAR-ARD/s1ard/compare/v1.3.0...v1.4.0>`_

1.3.0 | 2023-05-24
------------------

* SNAP RTC: increase DEM oversampling by a factor of two (`#78 <https://github.com/SAR-ARD/s1ard/pull/78>`_)
* nrb.format: do not hardcode src_nodata and read it from the data instead (`#79 <https://github.com/SAR-ARD/s1ard/pull/79>`_)
* enable configuration via command line arguments (`#80 <https://github.com/SAR-ARD/s1ard/pull/80>`_)
* improved date parsing (`#81 <https://github.com/SAR-ARD/s1ard/pull/81>`_)
* scene search via STAC (`#82 <https://github.com/SAR-ARD/s1ard/pull/82>`_)
* enhanced time filtering (`#84 <https://github.com/SAR-ARD/s1ard/pull/84>`_)
* general processor improvements (`#85 <https://github.com/SAR-ARD/s1ard/pull/85>`_)

`Full v1.3.0 Changelog <https://github.com/SAR-ARD/s1ard/compare/v1.2.0...v1.3.0>`_

1.2.0 | 2022-12-29
------------------

* improved geometry handling (`#71 <https://github.com/SAR-ARD/s1ard/pull/71>`_)
* DEM handling improvements (`#72 <https://github.com/SAR-ARD/s1ard/pull/72>`_)
* GRD buffering by (`#73 <https://github.com/SAR-ARD/s1ard/pull/73>`_)
* add DEM as additional output layer (`#70 <https://github.com/SAR-ARD/s1ard/pull/70>`_)
* sigma0 processing and annotation layer configuration (`#74 <https://github.com/SAR-ARD/s1ard/pull/74>`_)

`Full v1.2.0 Changelog <https://github.com/SAR-ARD/s1ard/compare/v1.1.0...v1.2.0>`_

1.1.0 | 2022-09-29
------------------

* documentation improvements (`#60 <https://github.com/SAR-ARD/s1ard/pull/60>`_)
* installation update (`#61 <https://github.com/SAR-ARD/s1ard/pull/61>`_)
* Process restructuring (`#63 <https://github.com/SAR-ARD/s1ard/pull/63>`_)
* minor structural changes and bug fixes (`#65 <https://github.com/SAR-ARD/s1ard/pull/65>`_)
* documentation update reflecting the recent process restructuring (`#66 <https://github.com/SAR-ARD/s1ard/pull/66>`_)
* renamed processing mode 'snap' to 'rtc' (`#67 <https://github.com/SAR-ARD/s1ard/pull/67>`_)

`Full v1.1.0 Changelog <https://github.com/SAR-ARD/s1ard/compare/v1.0.2...v1.1.0>`_

1.0.2 | 2022-08-24
------------------

* Fix error in handling of temporary VRTs (`#50 <https://github.com/SAR-ARD/s1ard/pull/50>`_)
* Adjustments to VRT log scaling (`#52 <https://github.com/SAR-ARD/s1ard/pull/52>`_)
* [metadata] read nodata values directly from files (instead of hard-coding them) (`#53 <https://github.com/SAR-ARD/s1ard/pull/53>`_)
* use type identifier in scene-specific DEM file names (`#55 <https://github.com/SAR-ARD/s1ard/pull/55>`_)
* Add VRT assets to STAC files (`#56 <https://github.com/SAR-ARD/s1ard/pull/56>`_)
* Fix and improve metadata geometry handling (`#57 <https://github.com/SAR-ARD/s1ard/pull/57>`_)
* SNAP 9 compatibility (`#58 <https://github.com/SAR-ARD/s1ard/pull/58>`_)

`Full v1.0.2 Changelog <https://github.com/SAR-ARD/s1ard/compare/v1.0.1...v1.0.2>`_

1.0.1 | 2022-07-03
------------------

* dem handling improvements (`#45 <https://github.com/SAR-ARD/s1ard/pull/45>`_)

`Full v1.0.1 Changelog <https://github.com/SAR-ARD/s1ard/compare/v1.0.0...v1.0.1>`_

1.0.0 | 2022-06-23
------------------

* Dockerfile to build s1ard image (`#31 <https://github.com/SAR-ARD/s1ard/pull/31>`_)
* adjustments to nodata value (`#28 <https://github.com/SAR-ARD/s1ard/pull/28>`_)
* renamed XML tag 'nrb' to 's1-nrb' (`#36 <https://github.com/SAR-ARD/s1ard/pull/36>`_)
* Metadata & Config Improvements (`#30 <https://github.com/SAR-ARD/s1ard/pull/30>`_)
* Geolocation accuracy (`#40 <https://github.com/SAR-ARD/s1ard/pull/40>`_)
* various bug fixes and documentation improvements

`Full v1.0.0 Changelog <https://github.com/SAR-ARD/s1ard/compare/v0.4.2...v1.0.0>`_

0.4.2 | 2022-06-16
------------------

* Update documentation (`#27 <https://github.com/SAR-ARD/s1ard/pull/27>`_)
* find unpacked .SAFE scenes in scene_dir (instead of just .zip) (`aea53a5 <https://github.com/SAR-ARD/s1ard/commit/aea53a57bc5fa1418fea4f46f69b41b7332909b1>`_)

`Full v0.4.2 Changelog <https://github.com/SAR-ARD/s1ard/compare/v0.4.1...v0.4.2>`_

0.4.1 | 2022-06-01
------------------

* handle ETAD products as zip, tar, and SAFE (`#25 <https://github.com/SAR-ARD/s1ard/pull/25>`_)
* set dem download authentication via env. variables (`#26 <https://github.com/SAR-ARD/s1ard/pull/26>`_)
* various bug fixes

`Full v0.4.1 Changelog <https://github.com/SAR-ARD/s1ard/compare/v0.4.0...v0.4.1>`_

0.4.0 | 2022-05-30
------------------

* outsourced and restructured DEM preparation functionality (`#18 <https://github.com/SAR-ARD/s1ard/pull/18>`_)
* outsourced ETAD correction to dedicated module (`#19 <https://github.com/SAR-ARD/s1ard/pull/19>`_)
* XML validation & improvements (`#17 <https://github.com/SAR-ARD/s1ard/pull/17>`_)
* Restructuring and cleanup (`#20 <https://github.com/SAR-ARD/s1ard/pull/20>`_)
* outsourced NRB formatting to dedicated module (`#21 <https://github.com/SAR-ARD/s1ard/pull/21>`_)
* extended acquisition mode support (`#22 <https://github.com/SAR-ARD/s1ard/pull/22>`_)
* Set up sphinx documentation (`#23 <https://github.com/SAR-ARD/s1ard/pull/23>`_)
* AOI scene selection (`#24 <https://github.com/SAR-ARD/s1ard/pull/24>`_)

`Full v0.4.0 Changelog <https://github.com/SAR-ARD/s1ard/compare/v0.3.0...v0.4.0>`_

0.3.0 | 2022-03-30
------------------

* Updated metadata module (`#9 <https://github.com/SAR-ARD/s1ard/pull/9>`_)
* Modified `prepare_dem` interface (`#10 <https://github.com/SAR-ARD/s1ard/pull/10>`_)
* Various improvements (`#11 <https://github.com/SAR-ARD/s1ard/pull/11>`_)
* Modified working directory structure (`#12 <https://github.com/SAR-ARD/s1ard/pull/12>`_)
* Updated `ancillary.py` (`#13 <https://github.com/SAR-ARD/s1ard/pull/13>`_)
* Added ETAD correction (`#14 <https://github.com/SAR-ARD/s1ard/pull/14>`_)
* Improved RGB composite (`#15 <https://github.com/SAR-ARD/s1ard/pull/15>`_)
* Store DEM/WBM tiles in UTM zones different to the native MGRS zone (`#16 <https://github.com/SAR-ARD/s1ard/pull/16>`_)

`Full v0.3.0 Changelog <https://github.com/SAR-ARD/s1ard/compare/v0.2.0...v0.3.0>`_

0.2.0 | 2022-03-03
------------------

`Full v0.2.0 Changelog <https://github.com/SAR-ARD/s1ard/compare/v0.1.0...v0.2.0>`_
