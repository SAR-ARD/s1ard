import os
import re
import time
import shutil
import tarfile as tf
from pyroSAR import identify
from spatialist.ancillary import finder
from s1etad_tools.cli.slc_correct import s1etad_slc_correct_main


def process(scene, etad_dir, out_dir, log):
    """
    Apply ETAD correction to a Sentinel-1 SLC product
    
    Parameters
    ----------
    scene: pyroSAR.ID
        the Sentinel-1 SLC scene
    etad_dir: str
        the directory containing ETAD products. This will be searched for products matching the defined SLC.
    out_dir: str
        the directory to store results. The ETAD product is unpacked to this directory if necessary.
        Two new sub-directories SLC_original SLC_ETAD and are created, which contain the original unpacked
        scene and the corrected one respectively.
    log: logging.Logger
        a logger object to write log info

    Returns
    -------
    pyroSAR.ID
        the corrected scene as pyroSAR ID object
    """
    slc_corrected_dir = os.path.join(out_dir, 'SLC_etad')
    os.makedirs(slc_corrected_dir, exist_ok=True)
    slc_base = os.path.basename(scene.scene).replace('.zip', '.SAFE')
    slc_corrected = os.path.join(slc_corrected_dir, slc_base)
    if not os.path.isdir(slc_corrected):
        start_time = time.time()
        acqtime = re.findall('[0-9T]{15}', os.path.basename(scene.scene))
        result = finder(etad_dir, ['_'.join(acqtime)], regex=True)
        try:
            if len(result) == 0:
                raise RuntimeError('cannot find ETAD product for scene {}'.format(scene.scene))
            
            if result[0].endswith('.tar'):
                etad_base = os.path.basename(result[0]).replace('.tar', '.SAFE')
                etad = os.path.join(out_dir, etad_base)
                if not os.path.isdir(etad):
                    archive = tf.open(result[0], 'r')
                    archive.extractall(out_dir)
                    archive.close()
            elif result[0].endswith('SAFE'):
                etad = result[0]
            else:
                raise RuntimeError('ETAD products are required to be .tar archives or .SAFE folders')
            scene.unpack(os.path.join(out_dir, 'SLC_original'), exist_ok=True)
            s1etad_slc_correct_main(s1_product=scene.scene,
                                    etad_product=etad,
                                    outdir=slc_corrected_dir,
                                    nthreads=2)
            shutil.rmtree(os.path.join(out_dir, 'SLC_original'))
            t = round((time.time() - start_time), 2)
            log.info('[   ETAD] -- {scene} -- {time}'.format(scene=scene.scene, time=t))
        except Exception as e:
            log.error('[   ETAD] -- {scene} -- {error}'.format(scene=scene.scene, error=e))
            raise
    else:
        msg = 'Already processed - Skip!'
        print('### ' + msg)
    return identify(slc_corrected)
