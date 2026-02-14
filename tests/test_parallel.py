"""Tests for hedgehog.utils.parallel utilities."""

from __future__ import annotations

from hedgehog.utils.parallel import parallel_map, resolve_n_jobs

# ---------------------------------------------------------------------------
# resolve_n_jobs
# ---------------------------------------------------------------------------


class TestResolveNJobs:
    """Tests for the resolve_n_jobs function."""

    def test_stage_config_takes_priority(self):
        assert (
            resolve_n_jobs(stage_config={"n_jobs": 4}, global_config={"n_jobs": 8}) == 4
        )

    def test_global_config_fallback(self):
        assert resolve_n_jobs(stage_config={}, global_config={"n_jobs": 8}) == 8

    def test_default_fallback(self):
        result = resolve_n_jobs(stage_config={}, global_config={}, default=2)
        assert result == 2

    def test_minus_one_uses_cpu_count(self, monkeypatch):
        monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)
        monkeypatch.setattr("hedgehog.utils.parallel.os.cpu_count", lambda: 16)
        assert resolve_n_jobs(default=-1) == 16

    def test_zero_uses_cpu_count(self, monkeypatch):
        monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)
        monkeypatch.setattr("hedgehog.utils.parallel.os.cpu_count", lambda: 8)
        assert resolve_n_jobs(stage_config={"n_jobs": 0}) == 8

    def test_slurm_env_var(self, monkeypatch):
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "32")
        assert resolve_n_jobs(default=-1) == 32

    def test_slurm_overrides_cpu_count(self, monkeypatch):
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "16")
        monkeypatch.setattr("hedgehog.utils.parallel.os.cpu_count", lambda: 64)
        assert resolve_n_jobs(default=-1) == 16

    def test_none_configs(self):
        # Should not raise, falls back to default
        result = resolve_n_jobs(stage_config=None, global_config=None, default=2)
        assert result == 2

    def test_minimum_is_one(self, monkeypatch):
        monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)
        monkeypatch.setattr("hedgehog.utils.parallel.os.cpu_count", lambda: None)
        assert resolve_n_jobs(default=-1) == 1


# ---------------------------------------------------------------------------
# parallel_map
# ---------------------------------------------------------------------------


def _square(x):
    """Top-level picklable function for testing."""
    return x * x


def _identity(x):
    return x


class TestParallelMap:
    """Tests for the parallel_map function."""

    def test_sequential_mode(self):
        result = parallel_map(_square, [1, 2, 3, 4], n_jobs=1)
        assert result == [1, 4, 9, 16]

    def test_parallel_mode(self):
        result = parallel_map(_square, list(range(20)), n_jobs=2)
        assert result == [x * x for x in range(20)]

    def test_empty_list(self):
        result = parallel_map(_square, [], n_jobs=4)
        assert result == []

    def test_single_item(self):
        result = parallel_map(_square, [5], n_jobs=4)
        assert result == [25]

    def test_results_match_sequential_and_parallel(self):
        items = list(range(100))
        sequential = parallel_map(_square, items, n_jobs=1)
        parallel = parallel_map(_square, items, n_jobs=2)
        assert sequential == parallel

    def test_custom_chunksize(self):
        result = parallel_map(_square, [1, 2, 3], n_jobs=2, chunksize=1)
        assert result == [1, 4, 9]

    def test_preserves_order(self):
        items = list(range(50))
        result = parallel_map(_identity, items, n_jobs=4)
        assert result == items
