import os
import re
import shutil
from math import ceil
import copy
from lxml import etree
import numpy as np
from pyproj import Geod
from osgeo import gdal, gdalconst
from scipy.interpolate import griddata
from datetime import datetime
from dateutil.parser import parse as dateparse
from spatialist.envi import HDRobject
from spatialist.ancillary import finder
from pyroSAR import identify, identify_many
from pyroSAR.ancillary import Lock, LockCollection
from pyroSAR.snap.auxil import (gpt, parse_recipe,
                                parse_node, orb_parametrize)
from s1ard.config import keyval_check
from cesard.tile_extraction import aoi_from_scene
from cesard.ancillary import datamask
from cesard.snap import geo, gsr, mli, rtc, sgr, postprocess

# do not remove, needed for interface
from cesard.snap import find_datasets, lsm_encoding, version_dict

import logging

log = logging.getLogger('s1ard')


# main interface

def config_to_string(config):
    """
    Convert the values of a configuration dictionary to strings.

    Parameters
    ----------
    config: dict
        the configuration as returned by :func:`get_config_section`

    Returns
    -------
    dict
        the dictionary with the same structure but values converted to strings.
    """
    out = {}
    allowed = get_config_keys()
    for k, v in config.items():
        if k == 'dem_prepare_mode':
            continue
        if k not in allowed:
            raise ValueError(f"key '{k}' not in allowed keys: {allowed}")
        if v is None:
            out[k] = 'None'
        elif k in ['allow_res_osv', 'clean_edges',
                   'clean_edges_pixels', 'cleanup']:
            out[k] = str(v)
        elif k == 'gpt_args' and isinstance(v, list):
            out[k] = ' '.join(v)
        else:
            out[k] = v
    return out


def get_config_keys():
    """
    Get all allowed configuration keys.
    
    Returns
    -------
    List[str]
    """
    return ['allow_res_osv', 'clean_edges', 'clean_edges_pixels', 'cleanup',
            'dem_resampling_method', 'gpt_args', 'img_resampling_method']


def get_config_section(parser, **kwargs):
    """
    Get the`config.ini` `SNAP` section content as a dictionary.
    
    Parameters
    ----------
    parser: configparser.ConfigParser
    kwargs: dict[str]

    Returns
    -------
    dict
    """
    out = {}
    defaults = {
        'allow_res_osv': 'True',
        'cleanup': 'True',
        'clean_edges': 'True',
        'clean_edges_pixels': '4',
        'dem_resampling_method': 'BILINEAR_INTERPOLATION',
        'gpt_args': 'None',
        'img_resampling_method': 'BILINEAR_INTERPOLATION',
    }
    if 'SNAP' in parser.sections():
        section = parser['SNAP']
        for k, v in defaults.items():
            if k not in section and k not in kwargs:
                section[k] = v
    else:
        parser.read_dict({'SNAP': defaults})
        section = parser['SNAP']
    
    kwargs_str = config_to_string(kwargs)
    for k, v in kwargs_str.items():
        if k in get_config_keys():
            section[k] = v
    
    for k, v in section.items():
        v = keyval_check(key=k, val=v, allowed_keys=get_config_keys())
        if k == 'gpt_args':
            if v is not None:
                v = v.split(' ')
        if k == 'clean_edges_pixels':
            v = section.getint(k)
        if k in ['allow_res_osv', 'clean_edges', 'cleanup']:
            v = section.getboolean(k)
        out[k] = v
    out['dem_prepare_mode'] = 'single-4326'
    return out


def get_metadata(scene, outdir):
    """
    Get processing metadata needed for ARD product metadata.

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
    dim = finder(scenedir, ['*.dim'])[0]
    scene = identify(dim)
    rlks, azlks = scene.meta['looks']
    return {'azlks': azlks,
            'rlks': rlks}


def process(scene, outdir, measurement, spacing, dem,
            dem_resampling_method='BILINEAR_INTERPOLATION',
            img_resampling_method='BILINEAR_INTERPOLATION',
            rlks=None, azlks=None, tmpdir=None, export_extra=None,
            allow_res_osv=True, clean_edges=True, clean_edges_pixels=4,
            neighbors=None, gpt_args=None, cleanup=True):
    """
    Main function for SAR processing with SNAP.

    Parameters
    ----------
    scene: str
        The SAR scene file name.
    outdir: str
        The output directory for storing the final results.
    measurement: {'sigma', 'gamma'}
        the backscatter measurement convention:

        - gamma: RTC gamma nought (:math:`\\gamma^0_T`)
        - sigma: RTC sigma nought (:math:`\\sigma^0_T`)
    spacing: int or float
        The output pixel spacing in meters.
    dem: str
        The DEM filename. Can be created with :func:`s1ard.dem.mosaic`.
    dem_resampling_method: str
        The DEM resampling method.
    img_resampling_method: str
        The image resampling method.
    rlks: int or None
        The number of range looks.
    azlks: int or None
        The number of azimuth looks.
    tmpdir: str or None
        Path to a temporary directory for intermediate products.
    export_extra: list[str] or None
        A list of ancillary layers to create. Default None: do not create any ancillary layers.
        Options:

         - DEM
         - gammaSigmaRatio: :math:`\\sigma^0_T / \\gamma^0_T`
         - sigmaGammaRatio: :math:`\\gamma^0_T / \\sigma^0_T`
         - incidenceAngleFromEllipsoid
         - layoverShadowMask
         - localIncidenceAngle
         - NESZ: noise equivalent sigma zero
         - projectedLocalIncidenceAngle
         - scatteringArea
         - lookDirection: range look direction angle
    allow_res_osv: bool
        Also allow the less accurate RES orbit files to be used?
    clean_edges: bool
        Erode noisy image edges? See :func:`pyroSAR.snap.auxil.erode_edges`.
        Does not apply to layover-shadow mask.
    clean_edges_pixels: int
        The number of pixels to erode.
    neighbors: list[str] or None
        (only applies to GRD) an optional list of neighboring scenes to add
        a buffer around the main scene using function :func:`grd_buffer`.
        If GRDs are processed compeletely independently, gaps are introduced
        due to a missing overlap. If `neighbors` is None or an empty list,
        buffering is skipped.
    gpt_args: list[str] or None
        a list of additional arguments to be passed to the gpt call

        - e.g. ``['-x', '-c', '2048M']`` for increased tile cache size and intermediate clearing
    cleanup: bool
        Delete intermediate files after successful process termination?

    Returns
    -------

    Examples
    --------
    >>> from s1ard import snap
    >>> scene = 'S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12.zip'
    >>> dem = 'S1A_IW_SLC__1SDV_20200103T170700_20200103T170727_030639_0382D5_6A12_DEM_EEA10.tif'
    >>> outdir = '.'
    >>> spacing = 10
    >>> rlks = 5
    >>> azlks = 1
    >>> export_extra = ['localIncidenceAngle', 'incidenceAngleFromEllipsoid',
    >>>                 'scatteringArea', 'layoverShadowMask', 'gammaSigmaRatio']
    >>> snap.process(scene=scene, outdir=outdir, spacing=spacing, dem=dem,
    >>>              rlks=rlks, azlks=azlks, export_extra=export_extra)
    """
    if measurement not in ['gamma', 'sigma']:
        raise RuntimeError("'measurement' must either be 'gamma' or 'sigma'")
    if export_extra is None:
        export_extra = []
    basename = os.path.splitext(os.path.basename(scene))[0]
    outdir_scene = os.path.join(outdir, basename)
    if tmpdir is None:
        tmpdir = outdir
        tmpdir_scene = os.path.join(tmpdir, basename + '_tmp')
    else:
        tmpdir_scene = os.path.join(tmpdir, basename)
    os.makedirs(outdir_scene, exist_ok=True)
    os.makedirs(tmpdir_scene, exist_ok=True)
    
    out_base = os.path.join(outdir_scene, basename)
    tmp_base = os.path.join(tmpdir_scene, basename)
    
    id = identify(scene)
    workflows = []
    
    apply_rtc = True
    
    apply_grd_buffering = neighbors is not None and len(neighbors) > 0
    ############################################################################
    # general pre-processing
    out_pre = tmp_base + '_pre.dim'
    out_pre_wf = out_pre.replace('.dim', '.xml')
    workflows.append(out_pre_wf)
    output_noise = 'NESZ' in export_extra
    with Lock(out_pre):
        if not os.path.isfile(out_pre):
            log.info('preprocessing main scene')
            pre(src=scene, dst=out_pre, workflow=out_pre_wf,
                allow_res_osv=allow_res_osv, output_noise=output_noise,
                output_beta0=apply_rtc, gpt_args=gpt_args,
                add_slice_num=apply_grd_buffering)
        else:
            log.info('main scene has already been preprocessed')
    ############################################################################
    # GRD buffering
    if apply_grd_buffering:
        # general preprocessing of neighboring scenes
        out_pre_neighbors = []
        for item in neighbors:
            basename_nb = os.path.splitext(os.path.basename(item))[0]
            tmpdir_nb = os.path.join(tmpdir, basename_nb)
            os.makedirs(tmpdir_nb, exist_ok=True)
            tmp_base_nb = os.path.join(tmpdir_nb, basename_nb)
            out_pre_nb = tmp_base_nb + '_pre.dim'
            out_pre_nb_wf = out_pre_nb.replace('.dim', '.xml')
            with Lock(out_pre_nb):
                if not os.path.isfile(out_pre_nb):
                    log.info(f'preprocessing GRD neighbor: {item}')
                    pre(src=item, dst=out_pre_nb, workflow=out_pre_nb_wf,
                        allow_res_osv=allow_res_osv, output_noise=output_noise,
                        output_beta0=apply_rtc, gpt_args=gpt_args)
                else:
                    log.info(f'GRD neighbor has already been preprocessed: {item}')
            out_pre_neighbors.append(out_pre_nb)
        ########################################################################
        # buffering
        out_buffer = tmp_base + '_buf.dim'
        out_buffer_wf = out_buffer.replace('.dim', '.xml')
        with Lock(out_buffer):
            if not os.path.isfile(out_buffer):
                log.info('buffering GRD scene with neighboring acquisitions')
                with LockCollection(out_pre_neighbors + [out_pre], soft=True):
                    try:
                        grd_buffer(src=out_pre, dst=out_buffer, workflow=out_buffer_wf,
                                   neighbors=out_pre_neighbors, gpt_args=gpt_args,
                                   buffer=10 * spacing)
                        workflows.append(out_buffer_wf)
                        out_pre = out_buffer
                    except RuntimeError:
                        log.info('did not perform buffering because the slice number '
                                 'could not be determined')
                        pass
            else:
                out_pre = out_buffer
                log.info('GRD scene has already been buffered')
    ############################################################################
    # range look direction angle
    if 'lookDirection' in export_extra:
        with Lock(out_pre):
            look_direction(dim=out_pre)
    ############################################################################
    # multi-looking
    out_mli = tmp_base + '_mli.dim'
    out_mli_wf = out_mli.replace('.dim', '.xml')
    with Lock(out_pre, soft=True):
        with Lock(out_mli):
            if not os.path.isfile(out_mli):
                mli(src=out_pre, dst=out_mli, workflow=out_mli_wf,
                    spacing=spacing, rlks=rlks, azlks=azlks, gpt_args=gpt_args)
    if not os.path.isfile(out_mli):
        out_mli = out_pre
    else:
        workflows.append(out_mli_wf)
    ############################################################################
    # radiometric terrain flattening
    out_rtc = out_gsr = out_sgr = None
    if apply_rtc:
        out_rtc = tmp_base + '_rtc.dim'
        out_rtc_wf = out_rtc.replace('.dim', '.xml')
        workflows.append(out_rtc_wf)
        output_sigma0_rtc = measurement == 'sigma' or 'gammaSigmaRatio' in export_extra
        with LockCollection([out_mli, dem], soft=True):
            with Lock(out_rtc):
                if not os.path.isfile(out_rtc):
                    log.info('radiometric terrain correction')
                    rtc(src=out_mli, dst=out_rtc, workflow=out_rtc_wf, dem=dem,
                        dem_resampling_method=dem_resampling_method,
                        sigma0=output_sigma0_rtc,
                        scattering_area='scatteringArea' in export_extra,
                        gpt_args=gpt_args)
        ########################################################################
        # gamma-sigma ratio computation
        out_gsr = None
        if 'gammaSigmaRatio' in export_extra:
            out_gsr = tmp_base + '_gsr.dim'
            out_gsr_wf = out_gsr.replace('.dim', '.xml')
            workflows.append(out_gsr_wf)
            with Lock(out_rtc, soft=True):
                with Lock(out_gsr):
                    if not os.path.isfile(out_gsr):
                        log.info('computing gamma-sigma ratio')
                        gsr(src=out_rtc, dst=out_gsr, workflow=out_gsr_wf,
                            gpt_args=gpt_args)
        ########################################################################
        # sigma-gamma ratio computation
        out_sgr = None
        if 'sigmaGammaRatio' in export_extra:
            out_sgr = tmp_base + '_sgr.dim'
            out_sgr_wf = out_sgr.replace('.dim', '.xml')
            workflows.append(out_sgr_wf)
            with Lock(out_rtc, soft=True):
                with Lock(out_sgr):
                    if not os.path.isfile(out_sgr):
                        log.info('computing sigma-gamma ratio')
                        sgr(src=out_rtc, dst=out_sgr, workflow=out_sgr_wf,
                            gpt_args=gpt_args)
    ############################################################################
    # geocoding
    
    # Process to multiple UTM zones or just one?
    # For testing purposes only.
    utm_multi = True
    
    def run():
        out_geo = out_base + '_geo_{}.dim'.format(epsg)
        out_geo_wf = out_geo.replace('.dim', '.xml')
        sources = list(filter(None, [out_mli, out_rtc, out_gsr, out_sgr]))
        with LockCollection(sources, soft=True):
            with Lock(out_geo):
                if not os.path.isfile(out_geo):
                    log.info(f'geocoding to EPSG:{epsg}')
                    scene1 = identify(out_mli)
                    pols = scene1.polarizations
                    bands0 = ['NESZ_{}'.format(pol) for pol in pols]
                    if measurement == 'gamma':
                        bands1 = ['Gamma0_{}'.format(pol) for pol in pols]
                    else:
                        bands0.extend(['Sigma0_{}'.format(pol) for pol in pols])
                        bands1 = []
                    if 'scatteringArea' in export_extra:
                        bands1.append('simulatedImage')
                    if 'lookDirection' in export_extra:
                        bands0.append('lookDirection')
                    geo(*sources,
                        dst=out_geo, workflow=out_geo_wf,
                        spacing=spacing, crs=epsg, geometry=ext,
                        export_extra=export_extra,
                        standard_grid_origin_x=align_x,
                        standard_grid_origin_y=align_y,
                        bands0=bands0, bands1=bands1, dem=dem,
                        dem_resampling_method=dem_resampling_method,
                        img_resampling_method=img_resampling_method,
                        gpt_args=gpt_args)
                    log.info('edge cleaning')
                    postprocess(out_geo, clean_edges=clean_edges,
                                clean_edges_pixels=clean_edges_pixels)
                    log.info('creating valid data masks')
                    out_geo_data = out_geo.replace('.dim', '.data')
                    pattern = r'(?:Gamma0|Sigma0)_[VH]{2}\.img$'
                    measurements = finder(out_geo_data, [pattern],
                                          regex=True, recursive=False)
                    dm_ras = os.path.join(out_geo_data, 'datamask.tif')
                    dm_vec = dm_ras.replace('.tif', '.gpkg')
                    dm_vec = datamask(measurement=measurements[0],
                                      dm_ras=dm_ras, dm_vec=dm_vec)
                else:
                    log.info(f'geocoding to EPSG:{epsg} has already been performed')
        for wf in workflows:
            wf_dst = os.path.join(outdir_scene, os.path.basename(wf))
            if wf != wf_dst and not os.path.isfile(wf_dst):
                shutil.copyfile(src=wf, dst=wf_dst)
    
    log.info('determining UTM zone overlaps')
    aois = aoi_from_scene(scene=id, multi=utm_multi)
    for aoi in aois:
        ext = aoi['extent']
        epsg = aoi['epsg']
        align_x = aoi['extent_utm']['xmin']
        align_y = aoi['extent_utm']['ymax']
        run()
    ############################################################################
    # delete intermediate files
    if cleanup:
        log.info('cleaning up')
        if id.product == 'GRD':
            # delete everything except *_pre.* products which are reused for buffering
            # this needs to be improved so that these products are also removed if they
            # are no longer needed for any buffering.
            items = finder(target=tmpdir_scene, matchlist=['*'],
                           foldermode=1, recursive=False)
            for item in items:
                if not re.search(r'_pre\.', item):
                    if os.path.isfile(item):
                        os.remove(item)
                    else:
                        shutil.rmtree(item)
        else:
            shutil.rmtree(tmpdir_scene)


def translate_annotation(annotation, measurement):
    """
    Translate s1ard annotation keys to SAR processor naming.
    
    Parameters
    ----------
    annotation: List[str]
        the s1ard annotation keys (e.g. ei, gs)
    measurement: str
        the SAR backscatter measurement (gamma|sigma)

    Returns
    -------
    List[str]
        the annotation layer keys as required by the SAR processor
    """
    export_extra = None
    lookup = {'dm': 'layoverShadowMask',
              'ei': 'incidenceAngleFromEllipsoid',
              'lc': 'scatteringArea',
              'ld': 'lookDirection',
              'li': 'localIncidenceAngle',
              'np': 'NESZ',
              'gs': 'gammaSigmaRatio',
              'sg': 'sigmaGammaRatio'}
    
    if annotation is not None:
        annotation = [
            'gs' if x == 'ratio' and measurement == 'gamma'
            else 'sg' if x == 'ratio'
            else x for x in annotation]
        export_extra = []
        for layer in annotation:
            # supported by the SAR processor
            if layer in lookup:
                export_extra.append(lookup[layer])
            # supported by the NRB processor
            elif layer in ['em', 'id']:
                continue
            else:
                log.warning(f'unsupported annotation layer: {layer}')
    return export_extra


# processor-specific functions

def grd_buffer(src, dst, workflow, neighbors, buffer=100, gpt_args=None):
    """
    GRD extent buffering.
    GRDs, unlike SLCs, do not overlap in azimuth.
    With this function, a GRD can be buffered using the neighboring acquisitions.
    First, all images are mosaicked using the `SliceAssembly` operator
    and then subsetted to the extent of the main scene including a buffer.
    The `SliceAssembly` operator needs info about the slice number (i.e., the
    ID/position inside the data take). If the value in the metadata is 0 (which
    can be the case in NRT slicing mode), the slice number is determined using
    function :func:`~s1ard.snap.nrt_slice_num`. If this fails, the function will
    raise an error.

    Parameters
    ----------
    src: str
        the file name of the source scene in BEAM-DIMAP format.
    dst: str
        the file name of the target scene. Format is BEAM-DIMAP.
    workflow: str
        the output SNAP XML workflow filename.
    neighbors: list[str]
        the file names of neighboring scenes
    buffer: int
        the buffer size in meters
    gpt_args: list[str] or None
        a list of additional arguments to be passed to the gpt call

        - e.g. ``['-x', '-c', '2048M']`` for increased tile cache size
          and intermediate clearing

    Raises
    ------
    RuntimeError
        if the slice number of a scene is 0, and it could not be determined
        from the acquisition time

    Returns
    -------

    """
    if len(neighbors) == 0:
        raise RuntimeError("the list of 'neighbors' is empty")
    
    scenes = identify_many([src] + neighbors, sortkey='start')
    wf = parse_recipe('blank')
    ############################################
    read_ids = []
    for scene in scenes:
        read = parse_node('Read')
        read.parameters['file'] = scene.scene
        wf.insert_node(read)
        read_ids.append(read.id)
    ############################################
    asm = parse_node('SliceAssembly')
    wf.insert_node(asm, before=read_ids)
    ############################################
    id_main = [x.scene for x in scenes].index(src)
    buffer_px = int(ceil(buffer / scenes[0].spacing[1]))
    xmin = 0
    width = scenes[id_main].samples
    if id_main == 0:
        ymin = 0
        height = scenes[id_main].lines + buffer_px
    else:
        ymin = scenes[0].lines - buffer_px
        factor = 1 if id_main == len(scenes) - 1 else 2
        height = scenes[id_main].lines + buffer_px * factor
    sub = parse_node('Subset')
    sub.parameters['region'] = [xmin, ymin, width, height]
    sub.parameters['geoRegion'] = ''
    sub.parameters['copyMetadata'] = True
    wf.insert_node(sub, before=asm.id)
    ############################################
    write = parse_node('Write')
    wf.insert_node(write, before=sub.id)
    write.parameters['file'] = dst
    write.parameters['formatName'] = 'BEAM-DIMAP'
    ############################################
    wf.write(workflow)
    gpt(xmlfile=workflow, tmpdir=os.path.dirname(dst),
        gpt_args=gpt_args)


def look_direction(dim):
    """
    Compute the per-pixel range look direction angle.
    This adds a new layer to an existing BEAM-DIMAP product.

    Steps performed:

    - read geolocation grid points
    - limit grid point list to those relevant to the image
    - for each point, compute the range direction angle to the next point in range direction.
    - interpolate the grid to the full image dimensions

    Notes
    -----
    - The interpolation depends on the location of the grid points relative to the image.
      Hence, by subsetting the image by an amount of pixels/lines different to the grid point
      sampling rate, the first and last points will no longer be in the first and last line respectively.
    - The list might get very large when merging the scene with neighboring acquisitions using
      SliceAssembly and this longer list significantly changes the interpolation result.
      The difference in interpolation can be mitigated by reducing the list of points to
      those inside the image and those just outside of it.

    Parameters
    ----------
    dim: str
        a BEAM-DIMAP metadata file (extension .dim)

    Returns
    -------

    """
    
    def interpolate(infile, method='linear'):
        with open(infile, 'rb') as f:
            tree = etree.fromstring(f.read())
        
        lats = []
        lons = []
        lines = []
        pixels = []
        rgtimes = []
        aztimes = []
        
        pols = tree.xpath("//MDElem[@name='standAloneProductInformation']"
                          "/MDATTR[@name='transmitterReceiverPolarisation']")
        polarization = pols[0].text.lower()
        re_ns = "http://exslt.org/regular-expressions"
        ann_pol = tree.xpath(f"//MDElem[@name='annotation']"
                             f"//MDElem[re:test(@name, '-{polarization}-', 'i')]",
                             namespaces={'re': re_ns})
        nlines = int(tree.find('Raster_Dimensions/NROWS').text)
        npixels = int(tree.find('Raster_Dimensions/NCOLS').text)
        abstract = tree.xpath("//MDElem[@name='Abstracted_Metadata']")[0]
        
        for ann in ann_pol:
            points = ann.xpath(".//MDElem[@name='geolocationGridPoint']")
            for point in points:
                pixel = int(point.find("./MDATTR[@name='pixel']").text)
                line = int(point.find("./MDATTR[@name='line']").text)
                lat = float(point.find("./MDATTR[@name='latitude']").text)
                lon = float(point.find("./MDATTR[@name='longitude']").text)
                rgtime = float(point.find("./MDATTR[@name='slantRangeTime']").text)
                aztime = dateparse(point.find("./MDATTR[@name='azimuthTime']").text)
                aztime = (aztime - datetime(1900, 1, 1)).total_seconds()
                pixels.append(pixel)
                lines.append(line)
                rgtimes.append(rgtime)
                aztimes.append(aztime)
                lats.append(lat)
                lons.append(lon)
        coords = list(zip(rgtimes, aztimes))
        
        flt = abstract.xpath("./MDATTR[@name='first_line_time']")[0].text
        flt = (dateparse(flt) - datetime(1900, 1, 1)).total_seconds()
        llt = abstract.xpath("./MDATTR[@name='last_line_time']")[0].text
        llt = (dateparse(llt) - datetime(1900, 1, 1)).total_seconds()
        lti = float(abstract.xpath("./MDATTR[@name='line_time_interval']")[0].text)
        
        # limit the coords to those relevant to the image
        # (SliceAssembly extends the list but a subsequent Subset does not shorten it)
        az_before = [x[1] for x in coords if x[1] < flt]
        tmp_min = (max(az_before) if len(az_before) > 0 else flt) - lti
        az_after = [x[1] for x in coords if x[1] > llt]
        tmp_max = (min(az_after) if len(az_after) > 0 else llt) + lti
        
        coords_sub = [x for x in coords if tmp_min <= x[1] <= tmp_max]
        
        values = []
        coords_select = []
        g = Geod(ellps='WGS84')
        for i, v in enumerate(coords):
            if v in coords_sub:
                if i + 1 < len(coords) and pixels[i] < pixels[i + 1]:
                    az12, az21, dist = g.inv(lons[i], lats[i], lons[i + 1], lats[i + 1])
                    values.append(az12)
                else:
                    az12, az21, dist = g.inv(lons[i], lats[i], lons[i - 1], lats[i - 1])
                    values.append(az21)
                coords_select.append(v)
        
        coords = np.array(coords_select)
        values = np.array(values)
        
        rgtime_fl_max = max([v for i, v in enumerate(rgtimes) if lines[i] == 0])
        rgtime_ll_max = max([v for i, v in enumerate(rgtimes) if lines[i] == max(lines)])
        rgtime_max = min([rgtime_fl_max, rgtime_ll_max])
        
        xi = np.linspace(min(rgtimes), rgtime_max, npixels)
        yi = np.linspace(flt, llt, nlines)
        xi, yi = np.meshgrid(xi, yi)
        zi = griddata(coords, values, (xi, yi), method=method, fill_value=0)
        return zi
    
    def write(array, out, reference, format='ENVI'):
        src_dataset = gdal.Open(reference)
        cols = src_dataset.RasterXSize
        rows = src_dataset.RasterYSize
        driver = gdal.GetDriverByName(format)
        
        dst_dataset = driver.Create(out, cols, rows, 1, gdalconst.GDT_Float32)
        
        dst_dataset.SetMetadata(src_dataset.GetMetadata())
        dst_dataset.SetGeoTransform(src_dataset.GetGeoTransform())
        dst_dataset.SetProjection(src_dataset.GetProjection())
        if format == 'GTiff':
            dst_dataset.SetGCPs(src_dataset.GetGCPs(), src_dataset.GetGCPSpatialRef())
        else:
            array = array.astype('float32').byteswap().newbyteorder()
        dst_band = dst_dataset.GetRasterBand(1)
        dst_band.WriteArray(array)
        dst_band.FlushCache()
        dst_band = None
        src_dataset = None
        dst_dataset = None
        driver = None
        if format == 'ENVI':
            hdrfile = os.path.splitext(out)[0] + '.hdr'
            with HDRobject(hdrfile) as hdr:
                hdr.byte_order = 1
                hdr.description = 'Sentinel-1 EW Level-1 GRD Product - Unit: deg'
                hdr.band_names = 'lookDirection'
                hdr.write()
            auxfile = out + '.aux.xml'
            os.remove(auxfile)
    
    def metadata(dim):
        with open(dim, 'rb') as f:
            root = etree.fromstring(f.read())
        
        image_interp = root.find('Image_Interpretation')
        ld_search = image_interp.xpath('Spectral_Band_Info[BAND_NAME="lookDirection"]')
        
        if len(ld_search) == 0:
            # number of bands
            element_bands = root.find('Raster_Dimensions/NBANDS')
            bands = int(element_bands.text)
            element_bands.text = str(bands + 1)
            
            # data access
            data_access = root.find('Data_Access')
            data_files = data_access.findall('./Data_File')
            last = data_files[-1]
            new = copy.deepcopy(last)
            fpath = new.find('DATA_FILE_PATH').attrib['href']
            fpath = fpath.replace(os.path.basename(fpath), 'lookDirection.hdr')
            new.find('DATA_FILE_PATH').attrib['href'] = fpath
            new.find('BAND_INDEX').text = str(bands)
            data_access.insert(data_access.index(last) + 1, new)
            
            # band info
            info = etree.SubElement(image_interp, 'Spectral_Band_Info')
            npixels = root.find('Raster_Dimensions/NCOLS').text
            nlines = root.find('Raster_Dimensions/NROWS').text
            etree.SubElement(info, 'BAND_INDEX').text = str(bands)
            etree.SubElement(info, 'BAND_DESCRIPTION').text = 'range look direction'
            etree.SubElement(info, 'BAND_NAME').text = 'lookDirection'
            etree.SubElement(info, 'BAND_RASTER_WIDTH').text = npixels
            etree.SubElement(info, 'BAND_RASTER_HEIGHT').text = nlines
            etree.SubElement(info, 'DATA_TYPE').text = 'float32'
            etree.SubElement(info, 'PHYSICAL_UNIT').text = 'deg'
            etree.SubElement(info, 'LOG10_SCALED').text = 'false'
            etree.SubElement(info, 'NO_DATA_VALUE_USED').text = 'true'
            etree.SubElement(info, 'NO_DATA_VALUE').text = '0.0'
            
            etree.indent(root, space='    ')
            tree = etree.ElementTree(root)
            tree.write(dim, pretty_print=True, xml_declaration=True, encoding='utf-8')
    
    data = dim.replace('.dim', '.data')
    out = os.path.join(data, 'lookDirection.img')
    if not os.path.isfile(out):
        log.info('look direction computation')
        ref = finder(target=data, matchlist=['*.img'])[0]
        arr = interpolate(infile=dim)
        write(array=arr, out=out, reference=ref)
        metadata(dim=dim)
    else:
        log.info('look direction has already been computed')


def nrt_slice_num(dim):
    """
    Check whether a product has a non-zero slice number and add it if not.
    In NRT Slicing mode, both `sliceNumber` and `totalSlices` are 0 in the manifest.safe file.
    `sliceNumber` is however needed in function :func:`~s1ard.snap.grd_buffer` for
    the SNAP operator `SliceAssembly`.
    The time from `segmentStartTime` to `last_line_time` is divided by
    the acquisition duration (`last_line_time` - `first_line_time`).
    `totalSlices` is set to 100, which is expected to exceed the maximum possible value.

    Parameters
    ----------
    dim: str
        the scene in BEAM-DIMAP format

    Raises
    ------
    RuntimeError
        if the slice number is 0, and it cannot be computed because
        the segment start time cannot be read from the metadata

    Returns
    -------

    """
    with open(dim, 'rb') as f:
        root = etree.fromstring(f.read())
    abstract = root.xpath("//MDElem[@name='Abstracted_Metadata']")[0]
    slice_num = abstract.xpath("./MDATTR[@name='slice_num']")[0]
    if slice_num.text == '0':
        flt = dateparse(abstract.xpath("./MDATTR[@name='first_line_time']")[0].text)
        llt = dateparse(abstract.xpath("./MDATTR[@name='last_line_time']")[0].text)
        try:
            sst = dateparse(root.xpath("//MDATTR[@name='segmentStartTime']")[0].text)
        except IndexError:
            raise RuntimeError('could not determine slice number '
                               'due to missing segment start time')
        aqd = llt - flt
        slice_num_new = str(int(round((llt - sst) / aqd)))
        slice_num.text = slice_num_new
        for item in root.xpath("//MDATTR[@name='sliceNumber']"):
            item.text = slice_num_new
        for item in root.xpath("//MDATTR[@name='totalSlices']"):
            item.text = '100'
        etree.indent(root, space='    ')
        tree = etree.ElementTree(root)
        tree.write(dim, pretty_print=True, xml_declaration=True,
                   encoding='utf-8')


def pre(src, dst, workflow, allow_res_osv=True, osv_continue_on_fail=False,
        output_noise=True, output_beta0=True, output_sigma0=True,
        output_gamma0=False, add_slice_num=True, gpt_args=None):
    """
    General SAR preprocessing. The following operators are used (optional steps in brackets):
    Apply-Orbit-File(->Remove-GRD-Border-Noise)->Calibration->ThermalNoiseRemoval(->TOPSAR-Deburst)

    Parameters
    ----------
    src: str
        the file name of the source scene
    dst: str
        the file name of the target scene. Format is BEAM-DIMAP.
    workflow: str
        the output SNAP XML workflow filename.
    allow_res_osv: bool
        Also allow the less accurate RES orbit files to be used?
    osv_continue_on_fail: bool
        Continue processing if no OSV file can be downloaded or raise an error?
    output_noise: bool
        output the noise power images?
    output_beta0: bool
        output beta nought backscatter needed for RTC?
    output_sigma0: bool
        output sigma nought backscatter needed for NESZ?
    output_gamma0: bool
        output gamma nought backscatter needed?
    add_slice_num: bool
        determine a slice number and add it to the product's metadata?
        This is only necessary if GRD buffering is intended.
        See :func:`~s1ard.snap.nrt_slice_num`.
    gpt_args: list[str] or None
        a list of additional arguments to be passed to the gpt call
        
        - e.g. ``['-x', '-c', '2048M']`` for increased tile cache size and intermediate clearing
    
    Returns
    -------

    See Also
    --------
    pyroSAR.snap.auxil.orb_parametrize
    """
    scene = identify(src)
    if not os.path.isfile(workflow):
        polarizations = scene.polarizations
        wf = parse_recipe('blank')
        ############################################
        read = parse_node('Read')
        read.parameters['file'] = scene.scene
        wf.insert_node(read)
        ############################################
        orb = orb_parametrize(scene=scene, formatName='SENTINEL-1',
                              allow_RES_OSV=allow_res_osv,
                              continueOnFail=osv_continue_on_fail)
        wf.insert_node(orb, before=read.id)
        last = orb
        ############################################
        if re.search('S1[A-Z]', scene.sensor) and scene.product == 'GRD':
            bn = parse_node('Remove-GRD-Border-Noise')
            wf.insert_node(bn, before=last.id)
            bn.parameters['selectedPolarisations'] = polarizations
            last = bn
        ############################################
        cal = parse_node('Calibration')
        wf.insert_node(cal, before=last.id)
        cal.parameters['selectedPolarisations'] = polarizations
        cal.parameters['outputBetaBand'] = output_beta0
        cal.parameters['outputSigmaBand'] = output_sigma0
        cal.parameters['outputGammaBand'] = output_gamma0
        ############################################
        tnr = parse_node('ThermalNoiseRemoval')
        wf.insert_node(tnr, before=cal.id)
        tnr.parameters['outputNoise'] = output_noise
        last = tnr
        ############################################
        if scene.product == 'SLC' and scene.acquisition_mode in ['EW', 'IW']:
            deb = parse_node('TOPSAR-Deburst')
            wf.insert_node(deb, before=last.id)
            last = deb
        ############################################
        write = parse_node('Write')
        wf.insert_node(write, before=last.id)
        write.parameters['file'] = dst
        write.parameters['formatName'] = 'BEAM-DIMAP'
        ############################################
        wf.write(workflow)
    if not os.path.isfile(dst):
        gpt(xmlfile=workflow, tmpdir=os.path.dirname(dst),
            gpt_args=gpt_args, removeS1BorderNoiseMethod='ESA')
        if scene.product == 'GRD' and add_slice_num:
            try:
                nrt_slice_num(dim=dst)
            except RuntimeError:
                raise RuntimeError('cannot obtain slice number')
