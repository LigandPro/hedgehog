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


_INIT_TOKEN = None


def _init_token(value: str) -> None:
    global _INIT_TOKEN
    _INIT_TOKEN = value


def _use_token(x: int) -> str:
    return f"{_INIT_TOKEN}:{x}"


@pytest.mark.skipif(
    os.name == "nt", reason="multiprocessing fork semantics differ on Windows"
)
def test_parallel_map_progress_multiprocessing_with_initializer():
    calls: list[tuple[int, int]] = []

    def progress(done: int, total: int) -> None:
        calls.append((done, total))

    out = parallel_map(
        _use_token,
        [1, 2, 3, 4],
        n_jobs=2,
        chunksize=1,
        progress=progress,
        initializer=_init_token,
        initargs=("t",),
    )
    assert out == ["t:1", "t:2", "t:3", "t:4"]
    assert len(calls) == 4
    assert calls[0] == (1, 4)
    assert calls[-1] == (4, 4)
