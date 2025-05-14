import re
import pytest
from pyroSAR.drivers import identify, Archive
from s1ard.search import ASFArchive, STACArchive, STACParquetArchive, collect_neighbors, scene_select
from shapely import wkt


@pytest.mark.parametrize('archive_class', [
    Archive, ASFArchive, STACArchive, STACParquetArchive
])
def test_archive(archive_class, stac, stac_parquet, testdata, tmpdir):
    """
    Test archive functionality for different archive types

    Parameters
    ----------
    archive_class : class
        The archive class to test (STACArchive or ASFArchive)
    stac : dict
        STAC configuration dictionary
    testdata : dict
        Test data dictionary
    """
    search_args = {'sensor': 'S1A', 'product': 'GRD',
                   'mindate': '20200708T182600', 'maxdate': '20200708T182800'}
    if archive_class is Archive:
        dbfile = str(tmpdir / 'scenes.db')
        archive_args = {'dbfile': dbfile, 'cleanup': False}
        file_ext = '.zip'
        with ASFArchive() as archive:
            scenes = archive.select(sensor='S1A', product='GRD',
                                    mindate='20200708T170000',
                                    maxdate='20200708T190000',
                                    return_value='ASF')
        with Archive(**archive_args) as archive:
            archive.insert(scenes)
    elif archive_class is ASFArchive:
        archive_args = {}
        file_ext = '.zip'
    elif archive_class is STACArchive:
        archive_args = {'url': stac['url'],
                        'collections': [stac['collection']]}
        search_args['check_exist'] = False
        file_ext = '.SAFE'
    elif archive_class is STACParquetArchive:
        archive_args = {'files': stac_parquet}
        file_ext = '.SAFE'
    else:
        raise RuntimeError(f'archive_class must be STACArchive or ASFArchive, '
                           f'is {type(archive_class)}')
    
    with archive_class(**archive_args) as archive:
        scenes = archive.select(**search_args)
        assert len(scenes) == 4
        
        # Test neighbor collection
        scene = identify(testdata['s1'])
        neighbors = collect_neighbors(archive=archive, scene=scene,
                                      stac_check_exist=False)
        pids = sorted([re.sub(f'{file_ext}$', '', x)[-4:]
                       for x in neighbors])
        assert pids == ['D160', 'DAAD']
        
        return_values = ['product', 'acquisition_mode', 'mindate', 'maxdate',
                         'sensor', 'frameNumber',
                         'geometry_wkt', 'geometry_wkb'
                         ]
        values = archive.select(**search_args, return_value=return_values)
        values = [dict(zip(return_values, x)) for x in values]
        
        assert values[0]['product'] == 'GRD'
        assert values[0]['acquisition_mode'] == 'IW'
        assert values[0]['mindate'] == '20200708T182614'
        assert values[0]['maxdate'] == '20200708T182643'
        assert values[0]['sensor'] == 'S1A'
        assert values[0]['frameNumber'] == '03DDAA'
        
        geom = wkt.loads('POLYGON (('
                         '-3.54075 4.290702, -1.313489 4.754753, '
                         '-1.666533 6.50591, -3.901151 6.046738, '
                         '-3.54075 4.290702))')
        assert wkt.loads(values[0]['geometry_wkt']) == geom
        assert values[0]['geometry_wkb'] == geom.wkb
        
        scenes, tiles = scene_select(archive=archive, **search_args)
        
        # four scenes matching the search result:
        # S1A_IW_GRDH_1SDV_20200708T182614_20200708T182643_033367_03DDAA_D160 (first scene of the data take)
        # S1A_IW_GRDH_1SDV_20200708T182643_20200708T182708_033367_03DDAA_9550
        # S1A_IW_GRDH_1SDV_20200708T182708_20200708T182733_033367_03DDAA_DAAD
        # S1A_IW_GRDH_1SDV_20200708T182733_20200708T182758_033367_03DDAA_888C
        
        # 31 MGRS tiles overlapping with the four scenes
        
        # one additional scene to fully cover the MGRS tiles:
        # the default is "date_strict=True" so initially this scene is not selected
        # because it exceeds the defined time range.
        # S1A_IW_GRDH_1SDV_20200708T182758_20200708T182823_033367_03DDAA_A793
        
        assert len(scenes) == 5
        assert len(tiles) == 31


def test_aoi_date():
    """
    make sure the combination of aoi_tiles and a date range does not return incomplete results
    """
    #
    with ASFArchive() as archive:
        scenes, tiles = scene_select(archive=archive,
                                     sensor='S1A', product='GRD',
                                     mindate='20200207T051836', maxdate='20200207T051902',
                                     aoi_tiles=['33TUM'])
    assert len(scenes) == 2
    assert len(tiles) == 1
