import os
import pytest
from s1ard.config import get_config, init, write
from s1ard.snap import get_config_keys as snap_get_keys


def test_config(tmpdir):
    with pytest.raises(ValueError):
        config = get_config()
    with pytest.raises(ValueError):
        config = get_config(work_dir=str(tmpdir), aoi_tiles='xyz')
    config = get_config(work_dir=str(tmpdir), db_file='scenes.db', aoi_tiles='32TNT')
    assert config['processing']['aoi_tiles'] == ['32TNT']
    assert config['processing']['db_file'] == os.path.join(str(tmpdir), 'scenes.db')


def test_config_snap(tmpdir):
    config = get_config(work_dir=str(tmpdir), db_file='scenes.db', aoi_tiles='32TNT')
    assert 'snap' in config.keys()
    snap_keys = list(config['snap'].keys())
    del snap_keys[snap_keys.index('dem_prepare_mode')]
    assert sorted(snap_keys) == snap_get_keys()
    expected = {
        'allow_res_osv': True,
        'cleanup': True,
        'clean_edges': True,
        'clean_edges_pixels': 4,
        'dem_resampling_method': 'BILINEAR_INTERPOLATION',
        'gpt_args': None,
        'img_resampling_method': 'BILINEAR_INTERPOLATION'
    }
    
    for key, value in expected.items():
        assert config['snap'][key] == value
    
    out = str(tmpdir / 'config.ini')
    write(config, out)
    config = get_config(config_file=out)
    assert 'snap' in config.keys()
    snap_keys = list(config['snap'].keys())
    del snap_keys[snap_keys.index('dem_prepare_mode')]
    assert sorted(snap_keys) == snap_get_keys()
    for key, value in expected.items():
        assert config['snap'][key] == value


def test_init(tmpdir):
    target = str(tmpdir / 'config.ini')
    # work_dir undefined
    with pytest.raises(ValueError):
        init(target=target)
    # no search option defined
    with pytest.raises(RuntimeError):
        init(target=target, work_dir=str(tmpdir))
    gpt_args = ['-J-Xmx100G', '-c', '75G', '-q', '30']
    init(target=target, work_dir=str(tmpdir), db_file='scenes.db',
         gpt_args=gpt_args)
    config = get_config(target)
    assert config['snap']['gpt_args'] == gpt_args
    # file already exists
    with pytest.raises(RuntimeError):
        init(target=target, work_dir=str(tmpdir), db_file='scenes.db')
