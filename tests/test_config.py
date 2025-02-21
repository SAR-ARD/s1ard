import os
import pytest
from s1ard.config import get_config

def test_config(tmpdir):
    with pytest.raises(ValueError):
        config = get_config()
    with pytest.raises(ValueError):
        config = get_config(work_dir=str(tmpdir), aoi_tiles='xyz')
    config = get_config(work_dir=str(tmpdir), db_file='scenes.db', aoi_tiles='32TNT')
    assert config['processing']['aoi_tiles'] == ['32TNT']
    assert config['processing']['db_file'] == os.path.join(str(tmpdir), 'scenes.db')
