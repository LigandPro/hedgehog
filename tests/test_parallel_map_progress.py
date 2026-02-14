import os

import pytest

from hedgehog.utils.parallel import parallel_map


def _inc(x: int) -> int:
    return x + 1


def test_parallel_map_progress_sequential():
    calls: list[tuple[int, int]] = []

    def progress(done: int, total: int) -> None:
        calls.append((done, total))

    out = parallel_map(_inc, [1, 2, 3], n_jobs=1, progress=progress)
    assert out == [2, 3, 4]
    assert calls == [(1, 3), (2, 3), (3, 3)]


@pytest.mark.skipif(
    os.name == "nt", reason="multiprocessing fork semantics differ on Windows"
)
def test_parallel_map_progress_multiprocessing():
    calls: list[tuple[int, int]] = []

    def progress(done: int, total: int) -> None:
        calls.append((done, total))

    out = parallel_map(_inc, [1, 2, 3, 4], n_jobs=2, chunksize=1, progress=progress)
    assert out == [2, 3, 4, 5]
    assert len(calls) == 4
    assert calls[0] == (1, 4)
    assert calls[-1] == (4, 4)
