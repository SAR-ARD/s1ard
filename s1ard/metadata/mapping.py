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

# 'z_error': Maximum error threshold on values for LERC* compression.
# Will be ignored if a compression algorithm is used that isn't related to LERC.
LERC_ERR_THRES = {
    'vv-g-lin': 1e-4,
    'vh-g-lin': 1e-4,
    'hh-g-lin': 1e-4,
    'hv-g-lin': 1e-4,
    'vv-s-lin': 1e-4,
    'vh-s-lin': 1e-4,
    'hh-s-lin': 1e-4,
    'hv-s-lin': 1e-4,
    'ei': 1e-3,
    'em': 1e-3,
    'dm': 0.0,
    'li': 1e-2,
    'lc': 0.1,
    'ld': 1e-3,
    'gs': 1e-4,
    'id': 0.0,
    'np-vv': 2e-5,
    'np-vh': 2e-5,
    'np-hh': 2e-5,
    'np-hv': 2e-5,
    'sg': 1e-4,
    'wm': 2e-5
}

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

DEM_MAP = {
    'GETASSE30':
        {'access': 'https://step.esa.int/auxdata/dem/GETASSE30',
         'ref': 'https://seadas.gsfc.nasa.gov/help-8.1.0/desktop/GETASSE30ElevationModel.html',
         'type': 'elevation',
         'gsd': '30 arcsec',
         'egm': 'https://apps.dtic.mil/sti/citations/ADA166519'},
    'Copernicus 10m EEA DEM':
        {'access': 'ftps://cdsdata.copernicus.eu/DEM-datasets/COP-DEM_EEA-10-DGED/2021_1',
         'ref': 'https://spacedata.copernicus.eu/web/cscda/dataset-details?articleId=394198',
         'type': 'surface',
         'gsd': '10 m',
         'egm': 'https://bgi.obs-mip.fr/data-products/grids-and-models/egm2008-global-model/'},
    'Copernicus 30m Global DEM':
        {'access': 'https://copernicus-dem-30m.s3.eu-central-1.amazonaws.com/',
         'ref': 'https://copernicus-dem-30m.s3.amazonaws.com/readme.html',
         'type': 'surface',
         'gsd': '30 m',
         'egm': 'https://bgi.obs-mip.fr/data-products/grids-and-models/egm2008-global-model/'},
    'Copernicus 30m Global DEM II':
        {'access': 'ftps://cdsdata.copernicus.eu/DEM-datasets/COP-DEM_GLO-30-DGED/2021_1',
         'ref': 'https://spacedata.copernicus.eu/web/cscda/dataset-details?articleId=394198',
         'type': 'surface',
         'gsd': '30 m',
         'egm': 'https://bgi.obs-mip.fr/data-products/grids-and-models/egm2008-global-model/'},
    'Copernicus 90m Global DEM II':
        {'access': 'ftps://cdsdata.copernicus.eu/DEM-datasets/COP-DEM_GLO-90-DGED/2021_1',
         'ref': 'https://spacedata.copernicus.eu/web/cscda/dataset-details?articleId=394198',
         'type': 'surface',
         'gsd': '90 m',
         'egm': 'https://bgi.obs-mip.fr/data-products/grids-and-models/egm2008-global-model/'}
}

# XML namespaces are identifiers, and it is not their goal to be directly usable for schema retrieval:
# https://stackoverflow.com/a/30761004
NS_MAP = {
    'placeholder': 'http://earth.esa.int/sentinel-1/spec/role/1.0',
    'sar': 'http://www.opengis.net/sar/2.1',
    'eop': 'http://www.opengis.net/eop/2.1',
    'om': 'http://www.opengis.net/om/2.0',
    'gml': 'http://www.opengis.net/gml/3.2',
    'ows': 'http://www.opengis.net/ows/2.0',
    'xlink': 'http://www.w3.org/1999/xlink'}

ASSET_MAP = {
    '-dm.tif': {'type': 'Mask',
                'unit': None,
                'role': 'data-mask',
                'title': 'Data Mask Image',
                'allowed': ['not layover, nor shadow',
                            'layover',
                            'shadow',
                            'layover and shadow',
                            'ocean',
                            'lakes',
                            'rivers']},
    '-ei.tif': {'type': 'Angle',
                'unit': 'deg',
                'role': 'ellipsoid-incidence-angle',
                'title': 'Ellipsoid Incidence Angle'},
    '-em.tif': {'type': 'Elevation',
                'unit': 'meters',
                'role': 'digital-elevation-model',
                'title': 'Digital Elevation Model'},
    '-lc.tif': {'type': 'Scattering Area',
                'unit': 'square_meters',
                'role': 'contributing-area',
                'title': 'Local Contributing Area'},
    '-ld.tif': {'type': 'Angle',
                'unit': 'deg',
                'role': 'range-look-direction-angle',
                'title': 'Range Look Direction Angle'},
    '-li.tif': {'type': 'Angle',
                'unit': 'deg',
                'role': 'local-incidence-angle',
                'title': 'Local Incidence Angle'},
    '-gs.tif': {'type': 'Ratio',
                'unit': None,
                'role': 'gamma-sigma-ratio',
                'title': 'Gamma0 RTC to sigma0 RTC ratio'},
    '-id.tif': {'type': 'AcqID',
                'unit': None,
                'role': 'acquisition-id',
                'title': 'Acquisition ID Image'},
    '-np-[vh]{2}.tif': {'type': 'Sigma-0',
                        'unit': 'dB',
                        'role': 'noise-power',
                        'title': 'Noise Power'},
    '-sg.tif': {'type': 'Ratio',
                'unit': None,
                'role': 'sigma-gamma-ratio',
                'title': 'Sigma0 RTC to gamma0 RTC ratio'},
    '-wm.tif': {'type': 'Sigma-0',
                'unit': None,
                'role': 'wind-modelled-backscatter',
                'title': 'wind-modelled backscatter (OCN CMOD NRCS)'}
}

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
    'ancillaryData_KML': 'https://sentinel.esa.int/documents/247904/1955685/S2A_OPER_GIP_TILPAR_MPC__'
                         '20151209T095117_V20150622T000000_21000101T000000_B00.kml',
    'card4l_nrb': 'https://ceos.org/ard/files/PFS/NRB/v5.5/CARD4L-PFS_NRB_v5.5.pdf',
    'card4l_orb': 'https://ceos.org/ard/files/PFS/ORB/v1.0/'
                  'CARD4L_Product_Family_Specification_Ocean_Radar_Backscatter-v1.0.pdf',
    'faradayRotationReference': None,
    'geoCorrAccuracyReference': f'https://s1ard.readthedocs.io/en/'
                                f'v{get_release_version(s1ard.__version__)}/general/geoaccuracy.html',
    'geoCorrAlgorithm': 'https://sentinel.esa.int/documents/247904/1653442/Guide-to-Sentinel-1-Geocoding.pdf',
    'griddingConventionURL': 'https://www.mgrs-data.org/data/documents/nga_mgrs_doc.pdf',
    'noiseRemovalAlgorithm': 'https://sentinel.esa.int/documents/247904/2142675/Thermal-Denoising-of-Products-Generated-by-Sentinel-1-IPF',
    'orbitDataAccess': 'https://step.esa.int/auxdata/orbits/Sentinel-1',
    'platformReference': {
        'sentinel-1a': 'https://database.eohandbook.com/database/missionsummary.aspx?missionID=575',
        'sentinel-1b': 'https://database.eohandbook.com/database/missionsummary.aspx?missionID=576',
        'sentinel-1c': 'https://database.eohandbook.com/database/missionsummary.aspx?missionID=577',
        'sentinel-1d': 'https://database.eohandbook.com/database/missionsummary.aspx?missionID=814'
    },
    'radiometricAccuracyReference': None,
    'RTCAlgorithm': 'https://doi.org/10.1109/Tgrs.2011.2120616',
    'sensorCalibration': 'https://sentinel.esa.int/web/sentinel/technical-guides/sentinel-1-sar/sar-instrument/calibration',
    'source_access': 'https://dataspace.copernicus.eu',
    'source_doi': 'https://sentinel.esa.int/documents/247904/1877131/Sentinel-1-Product-Specification',
    'windNormReferenceModel': "https://www.ecmwf.int/sites/default/files/3.1.pdf"
}
