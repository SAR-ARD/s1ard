import os
import re
import itertools
import shutil
from spatialist import bbox, Raster
from spatialist.envi import HDRobject
from spatialist.ancillary import finder
from pyroSAR import identify, identify_many
from pyroSAR.snap.auxil import gpt, parse_recipe, parse_node, \
    orb_parametrize, mli_parametrize, tc_parametrize, \
    sub_parametrize, erode_edges
from S1_NRB.tile_extraction import tiles_from_aoi, extract_tile
from S1_NRB.ancillary import get_max_ext


def mli(src, dst, workflow, spacing=None, rlks=None, azlks=None, allow_res_osv=True):
    """
    | create a multi-looked image (MLI) using the following operators (optional steps in brackets):
    | Apply-Orbit-File(->Remove-GRD-Border-Noise)->Calibration->ThermalNoiseRemoval(->TOPSAR-Deburst->Multilook)
    
    Parameters
    ----------
    src: str
        the file name of the source scene
    dst: str
        the file name of the target scene. Format is BEAM-DIMAP.
    workflow: str
        the output SNAP XML workflow filename.
    spacing: int or float
        the target pixel spacing for automatic determination of looks
        using function :func:`pyroSAR.ancillary.multilook_factors`.
        Overridden by arguments `rlks` and `azlks` if they are not None.
    rlks: int
        the number of range looks.
    azlks: int
        the number of azimuth looks.
    allow_res_osv: bool
        Also allow the less accurate RES orbit files to be used?

    Returns
    -------

    """
    scene = identify(src)
    polarizations = scene.polarizations
    wf = parse_recipe('blank')
    ############################################
    read = parse_node('Read')
    read.parameters['file'] = scene.scene
    wf.insert_node(read)
    ############################################
    orb = orb_parametrize(scene=scene, formatName='SENTINEL-1',
                          allow_RES_OSV=allow_res_osv)
    wf.insert_node(orb, before=read.id)
    last = orb
    ############################################
    if scene.sensor in ['S1A', 'S1B'] and scene.product == 'GRD':
        bn = parse_node('Remove-GRD-Border-Noise')
        wf.insert_node(bn, before=last.id)
        bn.parameters['selectedPolarisations'] = polarizations
        last = bn
    ############################################
    cal = parse_node('Calibration')
    wf.insert_node(cal, before=last.id)
    cal.parameters['selectedPolarisations'] = polarizations
    cal.parameters['outputBetaBand'] = True
    cal.parameters['outputSigmaBand'] = True
    cal.parameters['outputGammaBand'] = False
    ############################################
    tnr = parse_node('ThermalNoiseRemoval')
    wf.insert_node(tnr, before=cal.id)
    tnr.parameters['outputNoise'] = True
    last = tnr
    ############################################
    if scene.product == 'SLC' and scene.acquisition_mode in ['EW', 'IW']:
        deb = parse_node('TOPSAR-Deburst')
        wf.insert_node(deb, before=last.id)
        last = deb
    ############################################
    ml = mli_parametrize(scene=scene, spacing=spacing, rlks=rlks, azlks=azlks)
    if ml is not None:
        wf.insert_node(ml, before=last.id)
        last = ml
    ############################################
    write = parse_node('Write')
    wf.insert_node(write, before=last.id)
    write.parameters['file'] = dst
    write.parameters['formatName'] = 'BEAM-DIMAP'
    ############################################
    wf.write(workflow)
    gpt(xmlfile=workflow, tmpdir=os.path.dirname(dst))


def rtc(src, dst, workflow, dem, dem_resampling_method='BILINEAR_INTERPOLATION', sigma0=True, scattering_area=True):
    """
    Radiometric Terrain Flattening.
    
    Parameters
    ----------
    src: str
        the file name of the source scene
    dst: str
        the file name of the target scene. Format is BEAM-DIMAP.
    workflow: str
        the output SNAP XML workflow filename.
    dem: str
        the input DEM file name.
    dem_resampling_method: str
        the DEM resampling method.
    sigma0: bool
        output sigma0 RTC backscatter?
    scattering_area:bool
        output scattering area image?

    Returns
    -------

    """
    scene = identify(src)
    wf = parse_recipe('blank')
    ############################################
    read = parse_node('Read')
    read.parameters['file'] = scene.scene
    wf.insert_node(read)
    ############################################
    tf = parse_node('Terrain-Flattening')
    polarizations = scene.polarizations
    bands = ['Beta0_{}'.format(pol) for pol in polarizations]
    wf.insert_node(tf, before=read.id)
    tf.parameters['sourceBands'] = bands
    if 'reGridMethod' in tf.parameters.keys():
        tf.parameters['reGridMethod'] = False
    tf.parameters['outputSigma0'] = sigma0
    tf.parameters['outputSimulatedImage'] = scattering_area
    tf.parameters['demName'] = 'External DEM'
    tf.parameters['externalDEMFile'] = dem
    tf.parameters['externalDEMApplyEGM'] = False
    with Raster(dem) as ras:
        tf.parameters['externalDEMNoDataValue'] = ras.nodata
    tf.parameters['demResamplingMethod'] = dem_resampling_method
    last = tf
    ############################################
    write = parse_node('Write')
    wf.insert_node(write, before=last.id)
    write.parameters['file'] = dst
    write.parameters['formatName'] = 'BEAM-DIMAP'
    ############################################
    wf.write(workflow)
    gpt(xmlfile=workflow, tmpdir=os.path.dirname(dst))


def gsr(src, dst, workflow, src_sigma=None):
    """
    gamma-sigma ratio computation for either ellipsoidal and RTC sigma nought.
    
    Parameters
    ----------
    src: str
        the file name of the source scene. Both gamma and sigma bands are expected unless `src_sigma` is defined.
    dst: str
        the file name of the target scene. Format is BEAM-DIMAP.
    workflow: str
        the output SNAP XML workflow filename.
    src_sigma: str or None
        the optional file name of a second source scene from which to read the sigma bands.

    Returns
    -------

    """
    scene = identify(src)
    pol = scene.polarizations[0]
    wf = parse_recipe('blank')
    ############################################
    read = parse_node('Read')
    read.parameters['file'] = scene.scene
    wf.insert_node(read)
    last = read
    ############################################
    if src_sigma is not None:
        read.parameters['sourceBands'] = f'Gamma_{pol}'
        read2 = parse_node('Read')
        read2.parameters['file'] = src_sigma
        read2.parameters['sourceBands'] = f'Sigma_{pol}'
        wf.insert_node(read2)
        ########################################
        merge = parse_node('BandMerge')
        wf.insert_node(merge, before=[read.id, read2.id])
        last = merge
    ############################################
    math = parse_node('BandMaths')
    wf.insert_node(math, before=last.id)
    ratio = 'gammaSigmaRatio'
    expression = f'Sigma0_{pol} / Gamma0_{pol}'
    
    math.parameters.clear_variables()
    exp = math.parameters['targetBands'][0]
    exp['name'] = ratio
    exp['type'] = 'float32'
    exp['expression'] = expression
    exp['noDataValue'] = 0.0
    ############################################
    write = parse_node('Write')
    wf.insert_node(write, before=math.id)
    write.parameters['file'] = dst
    write.parameters['formatName'] = 'BEAM-DIMAP'
    ############################################
    wf.write(workflow)
    gpt(xmlfile=workflow, tmpdir=os.path.dirname(dst))


def geo(*src, dst, workflow, spacing, crs, geometry=None, buffer=0.01,
        export_extra=None, standard_grid_origin_x=0, standard_grid_origin_y=0,
        dem, dem_resampling_method='BILINEAR_INTERPOLATION',
        img_resampling_method='BILINEAR_INTERPOLATION', **bands):
    """
    geocoding
    
    Parameters
    ----------
    src
        variable number of input scene file names
    dst: str
        the file name of the target scene. Format is BEAM-DIMAP.
    workflow: str
        the target XML workflow file name
    spacing: int or float
        the target pixel spacing in meters
    crs: int or str
        the target coordinate reference system
    geometry: dict or spatialist.vector.Vector or str or None
        a vector geometry to limit the target product's extent
    buffer: int or float
        an additional buffer in degrees to add around `geometry`
    export_extra: list[str] or None
        a list of ancillary layers to write. Supported options:
        
         - DEM
         - incidenceAngleFromEllipsoid
         - layoverShadowMask
         - localIncidenceAngle
         - projectedLocalIncidenceAngle
    standard_grid_origin_x: int or float
        the X coordinate for pixel alignment
    standard_grid_origin_y: int or float
        the Y coordinate for pixel alignment
    dem: str
        the DEM file
    dem_resampling_method: str
        the DEM resampling method
    img_resampling_method: str
        the SAR image resampling method
    bands
        band ids for the input scenes in `args` as lists with keys bands<index>,
        e.g., bands1=['NESZ_VV'], bands2=['Gamma0_VV'], ...

    Returns
    -------

    """
    wf = parse_recipe('blank')
    ############################################
    scenes = identify_many(list(filter(None, src)))
    read_ids = []
    for i, scene in enumerate(scenes):
        read = parse_node('Read')
        read.parameters['file'] = scene.scene
        if f'bands{i}' in bands.keys():
            read.parameters['useAdvancedOptions'] = True
            read.parameters['sourceBands'] = bands[f'bands{i}']
        wf.insert_node(read)
        read_ids.append(read.id)
    ############################################
    if len(scenes) > 1:
        merge = parse_node('BandMerge')
        wf.insert_node(merge, before=read_ids)
        last = merge
    else:
        last = wf['Read']
    ############################################
    if geometry is not None:
        sub = sub_parametrize(scene=scenes[0], geometry=geometry, buffer=buffer)
        wf.insert_node(sub, before=last.id)
        last = sub
    ############################################
    tc = tc_parametrize(spacing=spacing, t_srs=crs,
                        export_extra=export_extra,
                        alignToStandardGrid=True,
                        externalDEMFile=dem,
                        externalDEMApplyEGM=False,
                        standardGridOriginX=standard_grid_origin_x,
                        standardGridOriginY=standard_grid_origin_y,
                        standardGridAreaOrPoint='area',
                        demResamplingMethod=dem_resampling_method,
                        imgResamplingMethod=img_resampling_method)
    wf.insert_node(tc, before=last.id)
    ############################################
    write = parse_node('Write')
    wf.insert_node(write, before=tc.id)
    write.parameters['file'] = dst
    write.parameters['formatName'] = 'BEAM-DIMAP'
    ############################################
    wf.write(workflow)
    gpt(xmlfile=workflow, tmpdir=os.path.dirname(dst))


def check_spacing(spacing):
    """
    perform a check whether the spacing fits into the MGRS tile boundaries
    
    Parameters
    ----------
    spacing: int or float
        the target pixel spacing in meters

    Returns
    -------

    """
    if 10980 % spacing != 0:
        raise RuntimeError(f'target spacing of {spacing} m does not align with MGRS tile size of 10980 m.')
    if 9780 % spacing != 0:
        raise RuntimeError(f'target spacing of {spacing} m does not align with MGRS tile overlap of 9780 m.')


def process(scene, outdir, spacing, kml, dem,
            dem_resampling_method='BILINEAR_INTERPOLATION',
            img_resampling_method='BILINEAR_INTERPOLATION',
            rlks=None, azlks=None, tmpdir=None, export_extra=None,
            allow_res_osv=True, slc_clean_edges=True, slc_clean_edges_pixels=4,
            cleanup=True):
    """
    
    Parameters
    ----------
    scene: str
    outdir: str
    spacing: int or float
    kml: str
    dem_type: str
    dem_dir: str or None
    dem_resampling_method: str
    img_resampling_method: str
    rlks: int or None
    azlks: int or None
    tmpdir: str or None
    export_extra: list[str] or None
    allow_res_osv: bool
    slc_clean_edges: bool
    slc_clean_edges_pixels: int
    cleanup: bool
    username: str or None
    password: str or None

    Returns
    -------

    """
    if tmpdir is None:
        tmpdir = outdir
    basename = os.path.splitext(os.path.basename(scene))[0]
    outdir_scene = os.path.join(outdir, basename)
    tmpdir_scene = os.path.join(tmpdir, basename)
    os.makedirs(outdir_scene, exist_ok=True)
    os.makedirs(tmpdir_scene, exist_ok=True)
    
    out_base = os.path.join(outdir_scene, basename)
    tmp_base = os.path.join(tmpdir_scene, basename)
    
    id = identify(scene)
    workflows = []
    ############################################################################
    # general pre-processing
    out_mli = tmp_base + '_mli.dim'
    out_mli_wf = out_mli.replace('.dim', '.xml')
    workflows.append(out_mli_wf)
    if not os.path.isfile(out_mli):
        mli(src=scene, dst=out_mli, workflow=out_mli_wf,
            spacing=spacing, rlks=rlks, azlks=azlks,
            allow_res_osv=allow_res_osv)
    ############################################################################
    # radiometric terrain flattening
    out_rtc = tmp_base + '_rtc.dim'
    out_rtc_wf = out_rtc.replace('.dim', '.xml')
    workflows.append(out_rtc_wf)
    if not os.path.isfile(out_rtc):
        rtc(src=out_mli, dst=out_rtc, workflow=out_rtc_wf, dem=dem,
            dem_resampling_method=dem_resampling_method,
            sigma0='gammaSigmaRatio' in export_extra,
            scattering_area='scatteringArea' in export_extra)
    ############################################################################
    # gamma-sigma ratio computation
    out_gsr = None
    if 'gammaSigmaRatio' in export_extra:
        out_gsr = tmp_base + '_gsr.dim'
        out_gsr_wf = out_gsr.replace('.dim', '.xml')
        workflows.append(out_gsr_wf)
        if not os.path.isfile(out_gsr):
            gsr(src=out_rtc, dst=out_gsr, workflow=out_gsr_wf)
    ############################################################################
    # geocoding
    with id.bbox() as geom:
        tiles = tiles_from_aoi(vectorobject=geom, kml=kml)
    
    for zone, group in itertools.groupby(tiles, lambda x: x[:2]):
        group = list(group)
        geometries = [extract_tile(kml=kml, tile=x) for x in group]
        epsg = geometries[0].getProjection(type='epsg')
        print(f'### processing EPSG:{epsg}')
        ext = get_max_ext(geometries=geometries)
        align_x = ext['xmin']
        align_y = ext['ymax']
        del geometries
        with bbox(coordinates=ext, crs=epsg) as geom:
            geom.reproject(projection=4326)
            ext = geom.extent
        out_geo = out_base + '_geo_{}.dim'.format(epsg)
        out_geo_wf = out_geo.replace('.dim', '.xml')
        if not os.path.isfile(out_geo):
            scene1 = identify(out_mli)
            pols = scene1.polarizations
            bands0 = ['NESZ_{}'.format(pol) for pol in pols]
            bands1 = ['Gamma0_{}'.format(pol) for pol in pols]
            if 'scatteringArea' in export_extra:
                bands1.append('simulatedImage')
            geo(out_mli, out_rtc, out_gsr, dst=out_geo, workflow=out_geo_wf,
                spacing=spacing, crs=epsg, geometry=ext,
                export_extra=export_extra,
                standard_grid_origin_x=align_x,
                standard_grid_origin_y=align_y,
                bands0=bands0, bands1=bands1, dem=dem,
                dem_resampling_method=dem_resampling_method,
                img_resampling_method=img_resampling_method)
            postprocess(out_geo, slc_clean_edges=slc_clean_edges,
                        slc_clean_edges_pixels=slc_clean_edges_pixels)
        for wf in workflows:
            wf_dst = os.path.join(outdir_scene, os.path.basename(wf))
            shutil.copyfile(src=wf, dst=wf_dst)
    if cleanup:
        shutil.rmtree(tmpdir_scene)


def postprocess(src, slc_clean_edges=True, slc_clean_edges_pixels=4):
    """
    performs SLC edge cleaning and sets the nodata value in the output ENVI HDR files.
    
    Parameters
    ----------
    src: str
        the file name of the source scene. Format is BEAM-DIMAP.
    slc_clean_edges: bool
        perform SLC edge cleaning?
    slc_clean_edges_pixels: int
        the number of pixels to erode during edge cleaning.

    Returns
    -------

    """
    if slc_clean_edges:
        erode_edges(src=src, only_boundary=True, pixels=slc_clean_edges_pixels)
    datadir = src.replace('.dim', '.data')
    hdrfiles = finder(target=datadir, matchlist=['*.hdr'])
    for hdrfile in hdrfiles:
        with HDRobject(hdrfile) as hdr:
            hdr.data_ignore_value = 0
            hdr.write(hdrfile)


def find_datasets(scene, outdir, epsg):
    """
    Find processed datasets for a scene in a certain CRS.
    
    Parameters
    ----------
    scene: str
        the file name of the SAR scene
    outdir: str
        the output directory in which to search for results
    epsg: int
        the EPSG code defining the output projection of the processed scenes.

    Returns
    -------
    dict or None
        Either None if no datasets were found or a dictionary with the
        following keys and values pointing to the file names
        (keys in brackets depending on available polarizations):
        
         - hh-g-lin: gamma nought RTC backscatter HH polarization
         - hv-g-lin: gamma nought RTC backscatter HV polarization
         - vh-g-lin: gamma nought RTC backscatter VH polarization
         - vv-g-lin: gamma nought RTC backscatter VV polarization
         - dm: layover-shadow data mask
         - ei: ellipsoidal incident angle
         - gs: gamma-sigma ratio
         - lc: local contributing area (aka scattering area)
         - li: local incident angle
         - (np-hh): noise power HH polarization
         - (np-hv): noise power HV polarization
         - (np-vh): noise power VH polarization
         - (np-vv): noise power VV polarization
    """
    basename = os.path.splitext(os.path.basename(scene))[0]
    scenedir = os.path.join(outdir, basename)
    subdir = os.path.join(scenedir, basename + f'_geo_{epsg}.data')
    if not os.path.isdir(subdir):
        raise RuntimeError('no RTC processing output found')
    lookup = {'dm': r'layoverShadowMask\.img',
              'ei': r'incidenceAngleFromEllipsoid\.img',
              'gs': r'gammaSigmaRatio_[VH]{2}\.img',
              'lc': r'simulatedImage_[VH]{2}\.img',
              'li': r'localIncidenceAngle\.img'}
    out = {}
    for key, pattern in lookup.items():
        match = finder(target=subdir, matchlist=[pattern], regex=True)
        if len(match) > 0:
            out[key] = match[0]
    pattern = r'Gamma0_(?P<pol>[VH]{2})\.img'
    gamma = finder(target=subdir, matchlist=[pattern], regex=True)
    for item in gamma:
        pol = re.search(pattern, item).group('pol')
        out[f'{pol.lower()}-g-lin'] = item
    pattern = r'NESZ_(?P<pol>[VH]{2})\.img'
    nesz = finder(target=subdir, matchlist=[pattern], regex=True)
    for item in nesz:
        pol = re.search(pattern, item).group('pol')
        out[f'np-{pol.lower()}'] = item
    if len(out) > 0:
        return out


def get_metadata(scene, outdir):
    """
    Get processing metadata needed for NRB metadata
    
    Parameters
    ----------
    scene: str
        the name of the SAR scene
    outdir: str
        the directory to search for processing output

    Returns
    -------
    dict
    """
    basename = os.path.splitext(os.path.basename(scene))[0]
    scenedir = os.path.join(outdir, basename)
    wf_mli = finder(scenedir, ['*mli.xml'])[0]
    wf = parse_recipe(wf_mli)
    if 'Multilook' in wf.operators:
        rlks = int(wf['Multilook'].parameters['nRgLooks'])
        azlks = int(wf['Multilook'].parameters['nAzLooks'])
    else:
        rlks = azlks = 1
    return {'azlks': azlks,
            'rlks': rlks}
