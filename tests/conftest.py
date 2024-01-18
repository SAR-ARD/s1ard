import os
import pytest
import requests


@pytest.fixture
def kml():
    remote = 'https://sentinel.esa.int/documents/247904/1955685/S2A_OPER_GIP_TILPAR_MPC__20151209T095117_V20150622T000000_21000101T000000_B00.kml'
    td = os.path.dirname(os.path.abspath(__file__))
    local_path = os.path.join(td, 'data')
    os.makedirs(local_path, exist_ok=True)
    local = os.path.join(local_path, os.path.basename(remote))
    if not os.path.isfile(local):
        print(f'downloading KML file to {local_path}')
        r = requests.get(remote)
        with open(local, 'wb') as out:
            out.write(r.content)
    return local


@pytest.fixture
def stac():
    return {'url': 'https://stac.terrabyte.lrz.de/public/api',
            'collection': 'sentinel-1-grd'}
