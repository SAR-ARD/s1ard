import pytest
import requests
from s1ard.metadata.mapping import URL


def url_recursive(key, dictionary, parent_key=None):
    key_info = f"{parent_key}.{key}" if parent_key else key
    if dictionary[key] is None:
        return
    if isinstance(dictionary[key], dict):
        for k, v in dictionary[key].items():
            url_recursive(k, dictionary[key], key_info)
    else:
        url = dictionary[key]
        try:
            response = requests.get(
                url=url,
                headers={"User-Agent": "Mozilla/5.0"},
                allow_redirects=True,
                timeout=30,
            )
        except requests.RequestException as exc:
            raise AssertionError(
                f'Could not reach {url}: {exc}'
            ) from exc
    
        assert response.status_code not in (404, 410), (
            f'URL no longer exists: {url}\n'
            f'Final URL: {response.url}\n'
            f'Status: {response.status_code}'
        )


@pytest.mark.parametrize('url_key', URL.keys())
def test_url(url_key, dictionary=URL):
    url_recursive(url_key, dictionary)
