from S1_NRB.search import STACArchive, ASFArchive, collect_neighbors, ASF


def test_stac(stac):
    archive = STACArchive(url=stac['url'],
                          collections=[stac['collection']])
    scenes = archive.select(sensor='S1A', mindate='20200708T182600',
                            maxdate='20200708T182800', check_exist=False)
    assert len(scenes) == 4


def test_asf():
    with ASFArchive() as archive:
        scenes = archive.select(sensor='S1A', product='GRD',
                                mindate='20200708T182600', maxdate='20200708T182800',
                                return_value=ASF)
        assert len(scenes) == 4
        neighbors = collect_neighbors(archive=archive, scene=scenes[0])
        assert len(neighbors) == 2
