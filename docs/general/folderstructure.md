# Folder Structure

The following describes the structure created to store intermediate and final files during a processor run.
The structure is based on the default configuration defined in the `config.ini` file and can be modified by a user.
Folders are highlighted in bold.

<summary><b>work_dir</b></summary>
<blockquote>

[//]: # (dem_dir)
<details>
<summary><b>dem_dir</b></summary>
<div class="admonition note">
<p class="admonition-title">Note</p>
DEM tiles in MGRS grid and WGS84 vertical datum for fast mosaicing of scene-specific DEMs during RTC processing.
Tiles with a non-native UTM zone additionally contain the EPSG code in the name.
For example, the native projection of tile 33TUL is 33N/EPSG:32633 but a variant in EPSG:32632 might exist for full coverage of a SAR scene.
</div>
<blockquote>

<details>
<summary><b>Copernicus 10m EEA DEM</b></summary>
<blockquote>

32TPR_DEM.tif  
32TPS_DEM.tif  
33TUL_32632_DEM.tif  
...

</blockquote>
</details>

<details>
<summary><b>Copernicus 30m Global DEM</b></summary>
<blockquote>

32TPR_DEM.tif  
32TPS_DEM.tif  
33TUL_32632_DEM.tif  
...

</blockquote>
</details>

<details>
<summary><b>Copernicus 30m Global DEM II</b></summary>
<blockquote>

32TPR_DEM.tif  
32TPS_DEM.tif  
33TUL_32632_DEM.tif  
...

</blockquote>
</details>

<details>
<summary><b>GETASSE30</b></summary>
<blockquote>

32TPR_DEM.tif  
32TPS_DEM.tif  
33TUL_32632_DEM.tif  
...

</blockquote>
</details>

</blockquote>
</details>

[//]: # (log_dir)
<details>
<summary><b>log_dir</b></summary>
<div class="admonition note">
<p class="admonition-title">Note</p>
Log files for each processor run containing the full processor configuration (<code class="docutils literal notranslate"><span class="pre">config.ini</span></code>), 
the versions of relevant installed software, and details on individual processing steps.
</div>
<blockquote>

20220719T1339_process.log  
20220719T1032_process.log  
20220708T1118_process.log  
...

</blockquote>
</details>

[//]: # (nrb_dir)
<details>
<summary><b>nrb_dir</b></summary>
<div class="admonition note">
<p class="admonition-title">Note</p>
The final S1-NRB tiles sorted into subfolders by MGRS tile.
Additional STAC files can be generated using function<a class="reference internal" href="../api.html#S1_NRB.metadata.stac.make_catalog" title="S1_NRB.metadata.stac.make_catalog"><code class="xref py py-func docutils literal notranslate"><span class="pre">S1_NRB.metadata.stac.make_catalog()</span></code></a>:
<ul>
  <li><code class="docutils literal notranslate"><span class="pre">collection.json</span></code>: a STAC collection file referencing the product-specific STAC item files per MGRS tile</li>
  <li><code class="docutils literal notranslate"><span class="pre">catalog.json</span></code>: a STAC catalog referencing all collections</li>
</ul> 
</div>
<blockquote>

<details>
<summary><b>32TPS</b></summary>
<blockquote>

[//]: # (NRB product)
<details>
<summary><b>S1A_IW_NRB__1SDV_20200103T170705_030639_0382D5_32TPS_8090</b></summary>
<blockquote>

[//]: # (annotation)
<details>
<summary><b>annotation</b></summary>
<blockquote>

s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-dm.tif  
s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-ei.tif  
s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-gs.tif  
s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-id.tif  
s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-lc.tif  
s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-li.tif  
s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-np-vh.tif  
s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-np-vv.tif

</blockquote>
</details>

[//]: # (measurement)
<details>
<summary><b>measurement</b></summary>
<blockquote>

s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-cc-g-lin.vrt  
s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-vh-g-lin.tif  
s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-vh-g-log.tif  
s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-vh-s-lin.tif  
s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-vh-s-log.tif  
s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-vv-g-lin.tif  
s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-vv-g-log.tif  
s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-vv-s-lin.tif  
s1a-iw-nrb-20200103t170705-030639-0382d5-32tps-vv-s-log.tif

</blockquote>
</details>

[//]: # (source)
<details>
<summary><b>source</b></summary>
<blockquote>

S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12.json
S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12.xml

</blockquote>
</details>

[//]: # (support)
<details>
<summary><b>support</b></summary>
<blockquote>

product.xsd  
source.xsd

</blockquote>
</details>

S1A_IW_NRB__1SDV_20200103T170705_030639_0382D5_32TPS_8090.json
S1A_IW_NRB__1SDV_20200103T170705_030639_0382D5_32TPS_8090.xml

</blockquote>
</details>

collection.json

</blockquote>
</details>

...  
catalog.json

</blockquote>
</details>

[//]: # (rtc_dir)
<details>
<summary><b>rtc_dir</b></summary>
<div class="admonition note">
<p class="admonition-title">Note</p>
The RTC processing output per source scene per UTM zone.
</div>
<blockquote>

[//]: # (SNAP output)
<details>
<summary><b>S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12</b></summary>
<blockquote>

<details>
<summary><b>32632</b></summary>
<blockquote>

S1A__IW___A_20200103T170700_Cal_NR_Deb_Orb_ML_TF_TC_proc.xml  
S1A__IW___A_20200103T170700_datamask.gpkg  
S1A__IW___A_20200103T170700_datamask.tif  
S1A__IW___A_20200103T170700_gammaSigmaRatio.tif  
S1A__IW___A_20200103T170700_incidenceAngleFromEllipsoid.tif  
S1A__IW___A_20200103T170700_layoverShadowMask.tif  
S1A__IW___A_20200103T170700_localIncidenceAngle.tif  
S1A__IW___A_20200103T170700_manifest.safe  
S1A__IW___A_20200103T170700_Orb_Cal_NR_Deb_ML_TC_proc.xml  
S1A__IW___A_20200103T170700_scatteringArea.tif  
S1A__IW___A_20200103T170700_VH_gamma0-rtc.tif  
S1A__IW___A_20200103T170700_VH_NESZ.tif  
S1A__IW___A_20200103T170700_VH_sigma0-rtc.tif  
S1A__IW___A_20200103T170700_VV_gamma0-rtc.tif  
S1A__IW___A_20200103T170700_VV_NESZ.tif  
S1A__IW___A_20200103T170700_VV_sigma0-rtc.tif  

</blockquote>
</details>
...

</blockquote>
</details>
...

</blockquote>
</details>

[//]: # (tmp_dir)
<details>
<summary><b>tmp_dir</b></summary>
<div class="admonition note">
<p class="admonition-title">Note</p>
Intermediate RTC processor files per scene per UTM zone. 
<ul>
  <li>EPSG code subfolder: scene-specific DEM mosaic and intermediate (SNAP) processor files</li>
  <li>unpacked ETAD files (*_ETA_*)</li>
  <li>SLC_etad subfolder: ETAD-corrected SLCs</li>
</ul> 

For example, the scene-specific DEM mosaic and intermediate SNAP files.
The latter are deleted after the process terminates.
</div>
<blockquote>

<details>
<summary><b>S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12</b></summary>
<blockquote>

<details>
<summary><b>32632</b></summary>
<blockquote>

S1A__IW___A_20200103T170700_EEA10_32632.tif

</blockquote>
</details>
...

<b>S1A_IW_ETA__AXDV_20200103T170700_20200103T170727_030639_0382D5_256B.SAFE</b>  
...

<details>
<summary><b>SLC_etad</b></summary>
<blockquote>

S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12.SAFE  
...

</blockquote>

</details>
</blockquote>
</details>
...
</blockquote>
</details>

[//]: # (wbm_dir)
<details>
<summary><b>wbm_dir</b></summary>
<div class="admonition note">
<p class="admonition-title">Note</p>
WBM tiles in MGRS grid and WGS84 vertical datum.
Tiles with a non-native UTM zone additionally contain the EPSG code in the name.
For example, The native projection of tile 33TUL is 33N/EPSG:32633 but a variant in EPSG:32632 might exist for full coverage of a SAR scene.
</div>
<blockquote>

<details>
<summary><b>Copernicus 10m EEA DEM</b></summary>
<blockquote>

32TPR_WBM.tif  
32TPS_WBM.tif  
33TUL_32632_WBM.tif  
...

</blockquote>
</details>

<details>
<summary><b>Copernicus 30m Global DEM II</b></summary>
<blockquote>

32TPR_WBM.tif  
32TPS_WBM.tif  
33TUL_32632_WBM.tif  
...

</blockquote>
</details>

</blockquote>
</details>

scenes.db

</blockquote>
