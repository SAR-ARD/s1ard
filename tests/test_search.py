from S1_NRB.search import STACArchive


def test_search(stac, kml):
    archive = STACArchive(url=stac['url'],
                          collections=[stac['collection']])
    scenes = archive.select(sensor='S1A', mindate='20200708T182600',
                            maxdate='20200708T182800', check_exist=False)
    assert len(scenes) == 4
