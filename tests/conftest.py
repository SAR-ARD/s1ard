import os
import pytest
import platform


@pytest.fixture(autouse=True)
def tmp_home(monkeypatch, tmp_path):
    home = tmp_path / 'tmp_home'
    home.mkdir()
    var = 'USERPROFILE' if platform.system() == 'Windows' else 'HOME'
    monkeypatch.setenv(var, str(home))
    assert os.path.expanduser('~') == str(home)
    yield home


@pytest.fixture
def stac():
    return {'url': 'https://stac.terrabyte.lrz.de/public/api',
            'collection': 'sentinel-1-grd'}
