import s1ard
from packaging.version import parse, Version

ARD_PATTERN = r'^(?P<sensor>S1[ABCD])_' \
              r'(?P<mode>IW|EW|S[1-6])_' \
              r'(?P<product>NRB|ORB)_' \
              r'(?P<resolution>_)' \
              r'(?P<processingLevel>1)' \
              r'(?P<category>S)' \
              r'(?P<pols>SH|SV|DH|DV)_' \
              r'(?P<start>[0-9]{8}T[0-9]{6})_' \
              r'(?P<orbitNumber>[0-9]{6})_' \
              r'(?P<dataTakeID>[0-9A-F]{6})_' \
              r'(?P<mgrsID>[0-9A-Z]{5})_' \
              r'(?P<ID>[0-9A-Z]{4})'

# Source data resolution
# https://sentinels.copernicus.eu/web/sentinel/technical-guides/sentinel-1-sar/products-algorithms/level-1-algorithms/single-look-complex
RES_MAP_SLC = {
    'IW': {'az': {'IW1': 22.5,
                  'IW2': 22.7,
                  'IW3': 22.6},
           'rg': {'IW1': 2.7,
                  'IW2': 3.1,
                  'IW3': 3.5}},
    'EW': {'az': {'EW1': 43.7,
                  'EW2': 44.3,
                  'EW3': 45.2,
                  'EW4': 45.6,
                  'EW5': 44.0},
           'rg': {'EW1': 7.9,
                  'EW2': 9.9,
                  'EW3': 11.6,
                  'EW4': 13.3,
                  'EW5': 14.4}},
    'SM': {'az': {'S1': 4.9,
                  'S2': 4.9,
                  'S3': 4.9,
                  'S4': 4.9,
                  'S5': 3.9,
                  'S6': 4.9},
           'rg': {'S1': 1.7,
                  'S2': 2.0,
                  'S3': 2.5,
                  'S4': 3.3,
                  'S5': 3.3,
                  'S6': 3.6}}
}

# https://sentinels.copernicus.eu/web/sentinel/technical-guides/sentinel-1-sar/products-algorithms/level-1-algorithms/ground-range-detected
_F = {
    'IW': None,
    'EW': None,
    'SM': {'az': {'S1': 8.1,
                  'S2': 8.7,
                  'S3': 8.8,
                  'S4': 8.9,
                  'S5': 8.9,
                  'S6': 9.1},
           'rg': {'S1': 8.1,
                  'S2': 8.4,
                  'S3': 8.8,
                  'S4': 9.0,
                  'S5': 9.2,
                  'S6': 9.2}}
}

_H = {
    'IW': {'az': {'IW1': 22.5,
                  'IW2': 22.6,
                  'IW3': 22.6},
           'rg': {'IW1': 20.4,
                  'IW2': 20.3,
                  'IW3': 20.5}},
    'EW': {'az': {'EW1': 51.5,
                  'EW2': 51.1,
                  'EW3': 51.3,
                  'EW4': 51.1,
                  'EW5': 51.5},
           'rg': {'EW1': 49.1,
                  'EW2': 50.3,
                  'EW3': 50.4,
                  'EW4': 50.7,
                  'EW5': 51.4}},
    'SM': {'az': {'S1': 21.3,
                  'S2': 23.0,
                  'S3': 23.2,
                  'S4': 23.6,
                  'S5': 23.4,
                  'S6': 24.1},
           'rg': {'S1': 21.4,
                  'S2': 22.2,
                  'S3': 23.1,
                  'S4': 23.7,
                  'S5': 24.2,
                  'S6': 24.4}}
}

_M = {
    'IW': {'az': {'IW1': 90.2,
                  'IW2': 90.6,
                  'IW3': 90.3},
           'rg': {'IW1': 87.9,
                  'IW2': 87.8,
                  'IW3': 88.7}},
    'EW': {'az': {'EW1': 90.1,
                  'EW2': 89.4,
                  'EW3': 86.9,
                  'EW4': 86.5,
                  'EW5': 90.1},
           'rg': {'EW1': 90.9,
                  'EW2': 93.1,
                  'EW3': 93.3,
                  'EW4': 93.8,
                  'EW5': 95.1}},
    'SM': {'az': {'S1': 77.3,
                  'S2': 83.8,
                  'S3': 84.4,
                  'S4': 85.7,
                  'S5': 85.1,
                  'S6': 87.5},
           'rg': {'S1': 78.0,
                  'S2': 80.8,
                  'S3': 84.1,
                  'S4': 86.4,
                  'S5': 87.9,
                  'S6': 88.8}}
}

RES_MAP_GRD = {
    'F': _F,
    'H': _H,
    'M': _M}

ENL_MAP_GRD = {
    'F': {'IW': None,
          'EW': None,
          'SM': (3.5, 3.5, 3.7, 3.5, 3.7, 3.5)},
    'H': {'IW': (4.4, 4.3, 4.3),
          'EW': (2.8, 2.7, 2.7, 2.7, 2.7),
          'SM': (26.8, 26.3, 29.7, 26.8, 29.7, 26.8)},
    'M': {'IW': (83.9, 81.2, 80.5),
          'EW': (15.2, 9.7, 9.6, 9.5, 9.6),
          'SM': (358.3, 350.5, 398.4, 358.3, 398.4, 358.3)}
}

OSV_MAP = {
    'PREORB': 'predicted',
    'RESORB': 'restituted',
    'POEORB': 'precise'}

# https://sentinel.esa.int/documents/247904/1653442/Guide-to-Sentinel-1-Geocoding.pdf
SLC_ACC_MAP = {
    'SM': {'ALE': {'rg': -3.02,
                   'az': 2.02},
           '1sigma': {'rg': 0.26,
                      'az': 0.41}},
    'IW': {'ALE': {'rg': -2.99,
                   'az': 2.03},
           '1sigma': {'rg': 0.22,
                      'az': 0.53}},
    'EW': {'ALE': {'rg': -3.41,
                   'az': 1.37},
           '1sigma': {'rg': 0.7,
                      'az': 2.27}}
}


def get_release_version(version_str: str) -> str:
    """
    Reduces a development version string to its last release version

    Parameters
    ----------
    version_str : str
        Version string (e.g. "v2.2.1.dev80+g34bba8c.d20250516")

    Returns
    -------
    str
        Release version (e.g. "2.2.0")
    """
    v = parse(version_str.lstrip('v'))
    if isinstance(v, Version):
        release = v.base_version
        parts = [int(x) for x in release.split('.')]
        if v.is_devrelease:
            if len(parts) >= 3:
                parts[-1] = max(0, parts[-1] - 1)
        return '.'.join(str(x) for x in parts)
    return version_str


URL = {
    'faradayRotationReference': None,
    'geoCorrAccuracyReference': f'https://s1ard.readthedocs.io/en/'
                                f'v{get_release_version(s1ard.__version__)}/general/geoaccuracy.html',
    'geoCorrAlgorithm': 'https://sentinel.esa.int/documents/247904/1653442/Guide-to-Sentinel-1-Geocoding.pdf',
    'noiseRemovalAlgorithm': 'https://sentinel.esa.int/documents/247904/2142675/Thermal-Denoising-of-Products-Generated-by-Sentinel-1-IPF',
    'orbitDataAccess': 'https://step.esa.int/auxdata/orbits/Sentinel-1',
    'radiometricAccuracyReference': None,
    'RTCAlgorithm': 'https://doi.org/10.1109/Tgrs.2011.2120616',
    'sensorCalibration': 'https://sentinel.esa.int/web/sentinel/technical-guides/sentinel-1-sar/sar-instrument/calibration',
    'source_access': 'https://dataspace.copernicus.eu',
    'source_doi': 'https://sentinel.esa.int/documents/247904/1877131/Sentinel-1-Product-Specification',
    'windNormReferenceModel': "https://www.ecmwf.int/sites/default/files/3.1.pdf"
}
