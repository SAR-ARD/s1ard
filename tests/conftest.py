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


@pytest.fixture
def testdir():
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')


@pytest.fixture
def testdata(testdir):
    out = {
        's1': os.path.join(testdir, 'S1A_IW_GRDH_1SDV_20150222T170750_20150222T170815_004739_005DD8_3768.zip'),
    }
    return out
