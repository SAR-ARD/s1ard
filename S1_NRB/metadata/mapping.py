NRB_PATTERN = r'^(?P<sensor>S1[AB])_' \
              r'(?P<mode>IW|EW)_' \
              r'(?P<product>NRB)_' \
              r'(?P<resolution>_)' \
              r'(?P<processingLevel>1)' \
              r'(?P<category>S)' \
              r'(?P<pols>SH|SV|DH|DV|VV|HH|HV|VH)_' \
              r'(?P<start>[0-9]{8}T[0-9]{6})_' \
              r'(?P<stop>[0-9]{8}T[0-9]{6})_' \
              r'(?P<orbitNumber>[0-9]{6})_' \
              r'(?P<dataTakeID>[0-9A-F]{6})_' \
              r'(?P<mgrsTile>[0-9A-Z]{5})'

SRC_PATTERN = r'^(?P<sensor>S1[AB])_' \
              r'(?P<mode>IW|EW|S[1-6]{1})_' \
              r'(?P<product>GRDH|SLC_)_' \
              r'(?P<processingLevel>1)' \
              r'(?P<category>S)' \
              r'(?P<pols>SH|SV|DH|DV|VV|HH|HV|VH)_' \
              r'(?P<start>[0-9]{8}T[0-9]{6})_' \
              r'(?P<stop>[0-9]{8}T[0-9]{6})_' \
              r'(?P<orbitNumber>[0-9]{6})_' \
              r'(?P<dataTakeID>[0-9A-F]{6})_' \
              r'(?P<ID>[0-9A-Z]{4})'

# Source data resolution
# https://sentinels.copernicus.eu/web/sentinel/technical-guides/sentinel-1-sar/products-algorithms/level-1-algorithms/single-look-complex
RES_MAP = {'IW': {'azimuthResolution': {'IW1': '22.5',
                                        'IW2': '22.7',
                                        'IW3': '22.6'},
                  'rangeResolution': {'IW1': '2.7',
                                      'IW2': '3.1',
                                      'IW3': '3.5'}},
           'EW': {'azimuthResolution': {'EW1': '43.7',
                                        'EW2': '44.3',
                                        'EW3': '45.2',
                                        'EW4': '45.6',
                                        'EW5': '44.0'},
                  'rangeResolution': {'EW1': '7.9',
                                      'EW2': '9.9',
                                      'EW3': '11.6',
                                      'EW4': '13.3',
                                      'EW5': '14.4'}},
           'SM': {'azimuthResolution': {'S1': '4.9',
                                        'S2': '4.9',
                                        'S3': '4.9',
                                        'S4': '4.9',
                                        'S5': '3.9',
                                        'S6': '4.9'},
                  'rangeResolution': {'S1': '1.7',
                                      'S2': '2.0',
                                      'S3': '2.5',
                                      'S4': '3.3',
                                      'S5': '3.3',
                                      'S6': '3.6'}}
           }

ORB_MAP = {'PREORB': 'predicted',
           'RESORB': 'restituted',
           'POEORB': 'precise'}

SAMPLE_MAP = {'-dm.tif': {'type': 'mask',
                          'unit': None,
                          'role': 'data-mask',
                          'title': 'Data Mask Image',
                          'values': {0: 'not layover, nor shadow',
                                     1: 'layover',
                                     2: 'shadow',
                                     3: 'layover and shadow',
                                     4: 'ocean water'}},
              '-ei.tif': {'type': 'angle',
                          'unit': 'deg',
                          'role': 'ellipsoid-incidence-angle',
                          'title': 'Ellipsoid Incidence Angle'},
              '-lc.tif': {'type': 'scattering area',
                          'unit': 'square_meters',
                          'role': 'contributing-area',
                          'title': 'Local Contributing Area'},
              '-li.tif': {'type': 'angle',
                          'unit': 'deg',
                          'role': 'local-incidence-angle',
                          'title': 'Local Incidence Angle'},
              '-gs.tif': {'type': 'ratio',
                          'unit': None,
                          'role': 'gamma-sigma-ratio',
                          'title': 'Gamma0 RTC to sigma0 RTC ratio'},
              '-id.tif': {'type': 'mask',
                          'unit': None,
                          'role': 'acquisition-id',
                          'title': 'Acquisition ID Image'},
              '-np-vv.tif': {'type': 'noise power VV',
                             'unit': None,
                             'role': 'noise-power',
                             'title': 'Noise Power VV'},
              '-np-vh.tif': {'type': 'noise power VH',
                             'unit': None,
                             'role': 'noise-power',
                             'title': 'Noise Power VH'},
              '-np-hh.tif': {'type': 'noise power HH',
                             'unit': None,
                             'role': 'noise-power',
                             'title': 'Noise Power HH'},
              '-np-hv.tif': {'type': 'noise power HV',
                             'unit': None,
                             'role': 'noise-power',
                             'title': 'Noise Power HV'}}
