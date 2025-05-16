import pytest
import requests
from s1ard.metadata.mapping import URL


def url_recursive(key, dictionary, parent_key=None):
    key_info = f"{parent_key}.{key}" if parent_key else key
    if isinstance(dictionary[key], dict):
        for k, v in dictionary[key].items():
            url_recursive(k, dictionary[key], key_info)
    else:
        response = requests.get(dictionary[key])
        assert response.status_code in [200, 418]


@pytest.mark.parametrize('url_key', URL.keys())
def test_url(url_key, dictionary=URL):
    url_recursive(url_key, dictionary)
