import os
import re
import time
import shutil
import tarfile as tf
import zipfile as zf
from pyroSAR import identify
from spatialist.ancillary import finder
from s1etad_tools.cli.slc_correct import s1etad_slc_correct_main
import logging

log = logging.getLogger('s1ard')


def process(scene, etad_dir, out_dir):
    """
    Apply ETAD correction to a Sentinel-1 SLC product.
    
    Parameters
    ----------
    scene: pyroSAR.drivers.ID
        The Sentinel-1 SLC scene.
    etad_dir: str
        The directory containing ETAD products. This will be searched for products matching the defined SLC.
    out_dir: str
        The directory to store results. The ETAD product is unpacked to this directory if necessary.
        Two new sub-directories SLC_original SLC_ETAD and are created, which contain the original unpacked
        scene and the corrected one respectively.

    Returns
    -------
    pyroSAR.drivers.ID
        The corrected scene as a pyroSAR ID object.
    """
    slc_corrected_dir = os.path.join(out_dir, 'SLC_etad')
    os.makedirs(slc_corrected_dir, exist_ok=True)
    slc_base = os.path.basename(scene.scene).replace('.zip', '.SAFE')
    slc_corrected = os.path.join(slc_corrected_dir, slc_base)
    if not os.path.isdir(slc_corrected):
        start_time = time.time()
        items = re.match(scene.pattern, os.path.basename(scene.file)).groupdict()
        pattern = '{sensor}_{beam}_ETA__AX{pols}_{start}_{stop}.*(SAFE|zip|tar)$'.format(**items)
        result = finder(etad_dir, [pattern], regex=True, foldermode=1)
        try:
            if len(result) == 0:
                raise RuntimeError('cannot find ETAD product for scene {}'.format(scene.scene))
            match = result[0]
            ext = os.path.splitext(match)[1]
            if ext in ['.tar', '.zip']:
                etad_base = os.path.basename(match).replace(ext, '.SAFE')
                etad = os.path.join(out_dir, etad_base)
                if not os.path.isdir(etad):
                    if ext == '.tar':
                        archive = tf.open(match, 'r')
                    else:
                        archive = zf.ZipFile(match, 'r')
                    archive.extractall(out_dir)
                    archive.close()
            elif ext == '.SAFE':
                etad = match
            else:
                raise RuntimeError('ETAD products are required to be .tar/.zip archives or .SAFE folders')
            scene.unpack(os.path.join(out_dir, 'SLC_original'), exist_ok=True)
            s1etad_slc_correct_main(s1_product=scene.scene,
                                    etad_product=etad,
                                    outdir=slc_corrected_dir,
                                    nthreads=2,
                                    order=0)  # using the default 1 introduces a bias of about -0.5 dB.
            shutil.rmtree(os.path.join(out_dir, 'SLC_original'))
            t = round((time.time() - start_time), 2)
            log.info(f'ETAD correction finished in {t} seconds')
        except Exception as e:
            log.error(e)
            raise
    else:
        msg = 'Already processed - Skip!'
        print('### ' + msg)
    return identify(slc_corrected)
