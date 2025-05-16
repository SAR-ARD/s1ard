import pytest
from s1ard.config import init


def test_init(tmpdir):
    target = str(tmpdir / 'config.ini')
    # work_dir undefined
    with pytest.raises(AssertionError):
        init(target=target)
    # no search option defined
    with pytest.raises(RuntimeError):
        init(target=target, work_dir=str(tmpdir))
    init(target=target, work_dir=str(tmpdir), db_file='scenes.db')
    # file already exists
    with pytest.raises(RuntimeError):
        init(target=target, work_dir=str(tmpdir), db_file='scenes.db')
