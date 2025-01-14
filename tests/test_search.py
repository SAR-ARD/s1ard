from pyroSAR import identify
from s1ard.search import STACArchive, ASFArchive, collect_neighbors, scene_select


def test_stac(stac, testdata):
    with STACArchive(url=stac['url'],
                     collections=[stac['collection']]) as archive:
        scenes = archive.select(sensor='S1A', mindate='20200708T182600',
                                maxdate='20200708T182800', check_exist=False)
        assert len(scenes) == 4
        scene = identify(testdata['s1'])
        neighbors = collect_neighbors(archive=archive, scene=scene,
                                      stac_check_exist=False)
        suffixes = sorted([x[-9:] for x in neighbors])
        assert suffixes == ['CEAB.SAFE', 'EBC3.SAFE']


def test_asf(testdata):
    with ASFArchive() as archive:
        scenes = archive.select(sensor='S1A', product='GRD',
                                mindate='20200708T182600', maxdate='20200708T182800',
                                return_value='ASF')
        assert len(scenes) == 4
        scene = identify(testdata['s1'])
        neighbors = collect_neighbors(archive=archive, scene=scene,
                                      stac_check_exist=False)
        suffixes = sorted([x[-8:] for x in neighbors])
        assert suffixes == ['CEAB.zip', 'EBC3.zip']


def test_scene_select():
    with ASFArchive() as archive:
        scenes, tiles = scene_select(archive=archive,
                                     sensor='S1A', product='GRD',
                                     mindate='20200708T182600', maxdate='20200708T182800')
    
    # four scenes matching the search result:
    # S1A_IW_GRDH_1SDV_20200708T182614_20200708T182643_033367_03DDAA_D160 (first scene of the data take)
    # S1A_IW_GRDH_1SDV_20200708T182643_20200708T182708_033367_03DDAA_9550
    # S1A_IW_GRDH_1SDV_20200708T182708_20200708T182733_033367_03DDAA_DAAD
    # S1A_IW_GRDH_1SDV_20200708T182733_20200708T182758_033367_03DDAA_888C
    
    # 31 MGRS tiles overlapping with the four scenes
    
    # one additional scene to fully cover the MGRS tiles:
    # the default is "date_strict=True" so initially this scene is not selected
    # because it exceeds the define time range.
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
