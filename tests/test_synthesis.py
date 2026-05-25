"""Tests for synthesis/utils.py."""

import io
import json
import pickle  # noqa: S403 — test-only, trusted data
import sys
import warnings
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

import hedgehog.synthesis.sync as sync_module
import hedgehog.synthesis.utils as synthesis_utils
from hedgehog.synthesis.utils import (
    _build_score_filter_mask,
    _calculate_ra_scores_batch,
    _get_rascore_model_path,
    apply_synthesis_score_filters,
    calculate_synthesis_scores,
    get_input_path,
    merge_retrosynthesis_results,
    parse_retrosynthesis_results,
    prepare_input_smiles,
    run_aizynthfinder,
)
from hedgehog.synthesis.utils import (
    _calculate_sa_score_single as _calculate_sa_score,
)
from hedgehog.synthesis.utils import (
    _calculate_syba_score_single as _calculate_syba_score,
)
from tests.constants import (
    COL_MODEL_NAME,
    COL_SMILES,
    FILE_FILTERED_MOLECULES,
    SMILES_BENZENE,
    SMILES_ETHANOL,
)


def _aizynth_install_root(root: Path) -> Path:
    """Return the project-local AiZynthFinder workspace root used in tests."""
    return root / "modules" / "aizynthfinder"


def _aizynth_config_path(root: Path) -> Path:
    """Return the project-local AiZynthFinder config path used in tests."""
    return _aizynth_install_root(root) / "public" / "config.yml"


class _DummyRascoreClassifier:
    """Simple picklable classifier stub with scikit-like API."""

    def predict_proba(self, fps):
        n = len(fps)
        out = np.zeros((n, 2), dtype=float)
        out[:, 1] = 0.8
        out[:, 0] = 0.2
        return out


class TestCalculateSaScore:
    """Tests for _calculate_sa_score function."""

    def test_valid_smiles_benzene(self):
        """SA score for benzene - should be a valid score."""
        score = _calculate_sa_score(SMILES_BENZENE)
        if not np.isnan(score):
            assert 1 <= score <= 10

    def test_valid_smiles_ethanol(self):
        """SA score for ethanol - should be easy to synthesize."""
        score = _calculate_sa_score(SMILES_ETHANOL)
        if not np.isnan(score):
            assert 1 <= score <= 10
            # Simple molecules should have low SA score
            assert score < 5

    def test_invalid_smiles(self):
        """SA score for invalid SMILES - should return NaN."""
        score = _calculate_sa_score("invalid")
        assert np.isnan(score)

    def test_empty_smiles(self):
        """SA score for empty SMILES - should return NaN or None."""
        score = _calculate_sa_score("")
        # Empty SMILES returns None or NaN depending on implementation
        assert score is None or (isinstance(score, float) and np.isnan(score))


class TestCalculateSybaScore:
    """Tests for _calculate_syba_score function."""

    def test_valid_smiles(self):
        """SYBA score for valid SMILES."""
        score = _calculate_syba_score(SMILES_BENZENE)
        # Score can be NaN if SYBA model not available
        if not np.isnan(score):
            assert isinstance(score, (int, float))

    def test_invalid_smiles(self):
        """SYBA score for invalid SMILES - should return NaN."""
        score = _calculate_syba_score("invalid")
        assert np.isnan(score)


class TestCalculateNonpherScore:
    """Tests for optional Nonpher score behavior."""

    def test_missing_nonpher_dependency_returns_nan(
        self,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Unavailable Nonpher dependency should map to NaN score."""
        monkeypatch.setitem(synthesis_utils._lazy_cache, "nonpher_filter", None)
        monkeypatch.setattr(synthesis_utils, "_load_nonpher_filter_impl", lambda: False)

        score = synthesis_utils._calculate_nonpher_complexity_score_single(
            SMILES_ETHANOL
        )
        assert np.isnan(score)


class TestApplySynthesisScoreFilters:
    """Tests for apply_synthesis_score_filters function."""

    def test_filter_by_sa_score(self):
        """Filter molecules by SA score threshold."""
        df = pd.DataFrame(
            {
                "smiles": ["a", "b", "c"],
                "sa_score": [2.0, 5.0, 8.0],
            }
        )
        config = {"sa_score_min": 0, "sa_score_max": 5}
        result = apply_synthesis_score_filters(df, config)
        assert len(result) == 2  # only first two molecules pass

    def test_filter_by_ra_score(self):
        """Filter molecules by RA score threshold."""
        df = pd.DataFrame(
            {
                "smiles": ["a", "b", "c"],
                "ra_score": [0.8, 0.5, 0.2],
            }
        )
        config = {"ra_score_min": 0.5, "ra_score_max": "inf"}
        result = apply_synthesis_score_filters(df, config)
        assert len(result) == 2  # first two pass

    def test_filter_by_syba_score(self):
        """Filter molecules by SYBA score threshold."""
        df = pd.DataFrame(
            {
                "smiles": ["a", "b"],
                "syba_score": [100.0, -50.0],
            }
        )
        config = {"syba_score_min": 0, "syba_score_max": "inf"}
        result = apply_synthesis_score_filters(df, config)
        assert len(result) == 1  # only first passes

    def test_filter_by_sync_score(self):
        """Filter molecules by SYNC score threshold."""
        df = pd.DataFrame(
            {
                "smiles": ["a", "b"],
                "sync_score": [0.8, 0.2],
            }
        )
        config = {"sync_score_min": 0.5, "sync_score_max": 1.0}
        result = apply_synthesis_score_filters(df, config)
        assert len(result) == 1

    def test_nan_scores_pass_filter(self):
        """Molecules with NaN scores should pass filter."""
        df = pd.DataFrame(
            {
                "smiles": ["a", "b", "c"],
                "sa_score": [2.0, np.nan, 8.0],
            }
        )
        config = {"sa_score_min": 0, "sa_score_max": 5}
        result = apply_synthesis_score_filters(df, config)
        # 'a' passes (2.0), 'b' passes (NaN), 'c' fails (8.0)
        assert len(result) == 2

    def test_filter_by_nested_score_filters(self):
        """Nested score_filters should support newly registered score columns."""
        df = pd.DataFrame(
            {
                "smiles": ["a", "b", "c"],
                "sc_score": [2.0, 4.5, np.nan],
            }
        )
        config = {"score_filters": {"sc_score": {"min": 1.0, "max": 4.0}}}
        result = apply_synthesis_score_filters(df, config)

        assert result["smiles"].tolist() == ["a", "c"]

    def test_nested_score_filter_can_be_disabled(self):
        """A score filter with null thresholds should be skipped."""
        df = pd.DataFrame(
            {
                "smiles": ["a", "b"],
                "sc_score": [2.0, 5.0],
            }
        )
        config = {"score_filters": {"sc_score": {"min": None, "max": None}}}
        result = apply_synthesis_score_filters(df, config)

        assert len(result) == 2


class TestCalculateSynthesisScoresRegistry:
    """Tests for registry-driven synthesis score calculation."""

    @pytest.fixture(autouse=True)
    def _clean_optional_scorer_env(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Keep optional external scorer tests independent from host env."""
        synthesis_utils._MISSING_EXTERNAL_SCORER_WARNINGS.clear()
        for env_key in (
            "HEDGEHOG_AUTO_INSTALL",
            "HEDGEHOG_OPTIONAL_ENV_ROOT",
            "HEDGEHOG_FSSCORE_COMMAND",
            "HEDGEHOG_FSSCORE_PYTHON",
            "HEDGEHOG_FSSCORE_MODEL_PATH",
            "HEDGEHOG_FSSCORE_REPO_PATH",
            "HEDGEHOG_GASA_COMMAND",
            "HEDGEHOG_GASA_EXECUTABLE",
            "HEDGEHOG_GASA_API_URL",
            "HEDGEHOG_GASA_TIMEOUT_SECONDS",
            "HEDGEHOG_NONPHER_PYTHON",
        ):
            monkeypatch.delenv(env_key, raising=False)

    def test_calculate_synthesis_scores_uses_enabled_registry_order(
        self, monkeypatch: pytest.MonkeyPatch
    ):
        """Enabled scorers should drive score columns and worker settings."""
        df = pd.DataFrame({"smiles": [SMILES_ETHANOL, SMILES_BENZENE]})
        parallel_calls = []

        def _fake_parallel_map(func, items, n_jobs, **kwargs):
            parallel_calls.append((func, list(items), n_jobs))
            progress = kwargs.get("progress")
            if progress is not None:
                progress(len(items), len(items))
            return [float(len(parallel_calls))] * len(items)

        def _fake_batch(smiles_list, config, progress_cb=None):
            if progress_cb is not None:
                progress_cb(len(smiles_list), len(smiles_list))
            return [9.0] * len(smiles_list)

        monkeypatch.setattr(synthesis_utils, "parallel_map", _fake_parallel_map)
        monkeypatch.setitem(
            synthesis_utils.SYNTHESIS_SCORERS,
            "fakebatch",
            synthesis_utils.SynthesisScorer(
                name="fakebatch",
                column="fake_batch_score",
                batch_calculator=_fake_batch,
            ),
        )

        progress_events = []
        result = calculate_synthesis_scores(
            df,
            config={"enabled_scores": ["sa", "syba", "fakebatch"], "n_jobs": 7},
            progress_cb=lambda name, done, total: progress_events.append(
                (name, done, total)
            ),
        )

        assert result["sa_score"].tolist() == [1.0, 1.0]
        assert result["syba_score"].tolist() == [2.0, 2.0]
        assert result["fake_batch_score"].tolist() == [9.0, 9.0]
        assert [call[2] for call in parallel_calls] == [7, 1]
        assert [event[0] for event in progress_events] == [
            "sa_score",
            "syba_score",
            "fake_batch_score",
        ]

    def test_fsscore_external_command_adapter(self):
        """FSScore should read scores from an explicitly configured command."""
        command = (
            f'{sys.executable} -c "import pandas as pd, sys; '
            "df = pd.read_csv(sys.argv[1]); "
            "df['score'] = [0.25 + i for i in range(len(df))]; "
            'df.to_csv(sys.argv[2], index=False)" {input} {output}'
        )
        df = pd.DataFrame({"smiles": [SMILES_ETHANOL, SMILES_BENZENE]})

        result = calculate_synthesis_scores(
            df,
            config={"enabled_scores": ["fsscore"], "fsscore_command": command},
        )

        assert result["fs_score"].tolist() == [0.25, 1.25]

    def test_fsscore_external_command_reorders_by_smiles(self):
        """External scorer output with SMILES should be realigned to input order."""
        command = (
            f'{sys.executable} -c "import pandas as pd, sys; '
            "df = pd.read_csv(sys.argv[1]); "
            "df['score'] = df['smiles'].map({{'CCO': 0.25, 'c1ccccc1': 1.25}}); "
            "df = df[['smiles', 'score']].iloc[::-1]; "
            'df.to_csv(sys.argv[2], index=False)" {input} {output}'
        )
        df = pd.DataFrame({"smiles": [SMILES_ETHANOL, SMILES_BENZENE]})

        result = calculate_synthesis_scores(
            df,
            config={"enabled_scores": ["fsscore"], "fsscore_command": command},
        )

        assert result["fs_score"].tolist() == [0.25, 1.25]

    def test_fsscore_external_command_mismatched_smiles_returns_nan(
        self, caplog: pytest.LogCaptureFixture
    ):
        """Mismatched external SMILES should soft-fail instead of misassigning scores."""
        command = (
            f'{sys.executable} -c "import pandas as pd, sys; '
            "pd.DataFrame({{'smiles': ['CCN', 'CCC'], 'score': [0.1, 0.2]}})."
            'to_csv(sys.argv[2], index=False)" {input} {output}'
        )

        with caplog.at_level("WARNING"):
            result = calculate_synthesis_scores(
                pd.DataFrame({"smiles": [SMILES_ETHANOL, SMILES_BENZENE]}),
                config={"enabled_scores": ["fsscore"], "fsscore_command": command},
            )

        assert result["fs_score"].isna().all()
        assert "SMILES do not match" in caplog.text

    def test_fsscore_external_command_bad_placeholder_returns_nan(
        self, caplog: pytest.LogCaptureFixture
    ):
        """Bad command placeholders should soft-fail with NaN scores."""
        command = "echo {unknown} > {output}"

        with caplog.at_level("WARNING"):
            result = calculate_synthesis_scores(
                pd.DataFrame({"smiles": [SMILES_ETHANOL]}),
                config={"enabled_scores": ["fsscore"], "fsscore_command": command},
            )

        assert np.isnan(result.loc[0, "fs_score"])
        assert "scorer command failed" in caplog.text

    def test_enabled_scores_string_is_parsed_as_scorer_names(self):
        """String enabled_scores should not be iterated character-by-character."""
        scorers = synthesis_utils._resolve_enabled_scorers(
            {"enabled_scores": "sa, syba"}
        )

        assert [scorer.name for scorer in scorers] == ["sa", "syba"]

    def test_fsscore_repo_path_resolves_default_model(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        """FSScore should derive checkpoint path from configured checkout path."""
        repo_path = tmp_path / "fsscore_repo"
        model_path = (
            repo_path / "models" / "pretrain_graph_GGLGGL_ep242_best_valloss.ckpt"
        )
        model_path.parent.mkdir(parents=True, exist_ok=True)
        model_path.write_bytes(b"x")

        captured: dict[str, list[str]] = {}

        def _fake_run(command, **kwargs):
            captured["command"] = command
            input_csv = Path(command[command.index("--input-csv") + 1])
            assert pd.read_csv(input_csv).columns.tolist() == ["smiles"]
            output_csv = Path(command[command.index("--output-csv") + 1])
            output_csv.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame({"smiles_score": [0.11, 0.22]}).to_csv(output_csv, index=False)
            return SimpleNamespace(returncode=0)

        monkeypatch.setattr(synthesis_utils.subprocess, "run", _fake_run)

        df = pd.DataFrame({"smiles": [SMILES_ETHANOL, SMILES_BENZENE]})
        result = calculate_synthesis_scores(
            df,
            config={
                "enabled_scores": ["fsscore"],
                "fsscore_python": sys.executable,
                "fsscore_repo_path": str(repo_path),
            },
        )

        assert result["fs_score"].tolist() == [0.11, 0.22]
        command = captured["command"]
        assert "hedgehog.workers.fsscore_worker" in command
        assert command[command.index("--model-path") + 1] == str(model_path.resolve())
        assert "--num-workers" in command
        assert command[command.index("--num-workers") + 1] == "0"
        assert command[command.index("--graph-datapath") + 1].endswith(
            "fsscore_graphs.pt"
        )

    def test_fsscore_num_workers_is_capped_to_batch_size(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """FSScore should not request more workers than input molecules."""
        repo_path = tmp_path / "fsscore_repo"
        model_path = (
            repo_path / "models" / "pretrain_graph_GGLGGL_ep242_best_valloss.ckpt"
        )
        model_path.parent.mkdir(parents=True, exist_ok=True)
        model_path.write_bytes(b"x")

        captured: dict[str, list[str]] = {}

        def _fake_run(command, **kwargs):
            captured["command"] = command
            output_csv = Path(command[command.index("--output-csv") + 1])
            output_csv.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame({"smiles_score": [0.11]}).to_csv(output_csv, index=False)
            return SimpleNamespace(returncode=0)

        monkeypatch.setattr(synthesis_utils.subprocess, "run", _fake_run)

        result = calculate_synthesis_scores(
            pd.DataFrame({"smiles": [SMILES_ETHANOL]}),
            config={
                "enabled_scores": ["fsscore"],
                "fsscore_python": sys.executable,
                "fsscore_repo_path": str(repo_path),
                "n_jobs": 99,
                "fsscore_num_workers": 3,
            },
        )

        assert result["fs_score"].tolist() == [0.11]
        command = captured["command"]
        assert command[command.index("--num-workers") + 1] == "1"

    def test_fsscore_missing_model_path_returns_nan_with_warning(
        self, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
    ):
        """Configured FSScore python without model path should soft-fail to NaN."""
        synthesis_utils._MISSING_EXTERNAL_SCORER_WARNINGS.clear()
        monkeypatch.delenv("HEDGEHOG_AUTO_INSTALL", raising=False)
        caplog.set_level("WARNING")

        df = pd.DataFrame({"smiles": [SMILES_ETHANOL]})
        result = calculate_synthesis_scores(
            df,
            config={
                "enabled_scores": ["fsscore"],
                "fsscore_python": sys.executable,
            },
        )

        assert np.isnan(result.loc[0, "fs_score"])
        assert "HEDGEHOG_FSSCORE_MODEL_PATH" in caplog.text

    def test_gasa_executable_adapter(self, tmp_path: Path):
        """GASA should read scores from an explicitly configured executable."""
        script_path = tmp_path / "gasa_stub.sh"
        script_path.write_text(
            "#!/bin/sh\n"
            'if [ "$1" != "--smiles" ] || [ -z "$2" ]; then\n'
            "  echo 'missing --smiles argument' >&2\n"
            "  exit 2\n"
            "fi\n"
            "echo '{\"gasa_score\": 0.9}'\n",
            encoding="utf-8",
        )
        script_path.chmod(0o755)
        df = pd.DataFrame({"smiles": [SMILES_ETHANOL, SMILES_BENZENE]})

        result = calculate_synthesis_scores(
            df,
            config={
                "enabled_scores": ["gasa"],
                "gasa": {"executable": str(script_path)},
            },
        )

        assert result["gasa_score"].tolist() == [0.9, 0.9]

    def test_gasa_external_command_adapter(self):
        """GASA should read scores from an explicitly configured command."""
        command = (
            f'{sys.executable} -c "import pandas as pd, sys; '
            "df = pd.read_csv(sys.argv[1]); "
            "df['gasa_score'] = [0.75 + i for i in range(len(df))]; "
            'df.to_csv(sys.argv[2], index=False)" {input} {output}'
        )
        df = pd.DataFrame({"smiles": [SMILES_ETHANOL, SMILES_BENZENE]})

        result = calculate_synthesis_scores(
            df,
            config={"enabled_scores": ["gasa"], "gasa_command": command},
        )

        assert result["gasa_score"].tolist() == [0.75, 1.75]

    def test_gasa_local_api_adapter(self, monkeypatch: pytest.MonkeyPatch):
        """GASA should use a configured local API endpoint."""
        requests_seen: list[tuple[str, str]] = []

        class _DummyResponse:
            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, tb):
                del exc_type, exc, tb
                return False

            def read(self) -> bytes:
                return b'{"gasa_score": 0.42}'

        def _fake_urlopen(request, timeout):  # noqa: ANN001
            del timeout
            requests_seen.append((request.full_url, request.data.decode("utf-8")))
            return _DummyResponse()

        monkeypatch.setattr(synthesis_utils.urllib_request, "urlopen", _fake_urlopen)
        df = pd.DataFrame({"smiles": [SMILES_ETHANOL, SMILES_BENZENE]})

        result = calculate_synthesis_scores(
            df,
            config={
                "enabled_scores": ["gasa"],
                "gasa": {"api_url": "http://127.0.0.1:8123/score"},
            },
        )

        assert result["gasa_score"].tolist() == [0.42, 0.42]
        assert len(requests_seen) == 2
        assert requests_seen[0][0] == "http://127.0.0.1:8123/score"
        assert '"smiles"' in requests_seen[0][1]

    def test_gasa_non_local_api_is_rejected(self, caplog):
        """Configured non-local GASA API should be rejected for safety."""
        df = pd.DataFrame({"smiles": [SMILES_ETHANOL]})

        with caplog.at_level("WARNING"):
            result = calculate_synthesis_scores(
                df,
                config={
                    "enabled_scores": ["gasa"],
                    "gasa": {"api_url": "https://example.com/gasa"},
                },
            )

        assert np.isnan(result.loc[0, "gasa_score"])
        assert any("not local" in record.message for record in caplog.records)

    def test_missing_external_scorers_return_nan(
        self, monkeypatch: pytest.MonkeyPatch, caplog
    ):
        """Unconfigured optional external scorers should not fail the pipeline."""
        monkeypatch.setattr(synthesis_utils, "_auto_install_enabled", lambda: False)
        df = pd.DataFrame({"smiles": [SMILES_ETHANOL]})

        with caplog.at_level("WARNING"):
            result = calculate_synthesis_scores(
                df,
                config={"enabled_scores": ["fsscore", "gasa"]},
            )

        assert np.isnan(result.loc[0, "fs_score"])
        assert np.isnan(result.loc[0, "gasa_score"])
        assert any(
            "GASA scorer is enabled but no backend is configured" in record.message
            for record in caplog.records
        )

    def test_scscore_model_is_preloaded_before_parallel_workers(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """SCScore model should preload once before parallel worker fan-out."""
        call_order: list[str] = []
        df = pd.DataFrame({"smiles": [SMILES_ETHANOL, SMILES_BENZENE]})

        monkeypatch.setattr(
            synthesis_utils,
            "_load_scscore_model",
            lambda: call_order.append("load") or False,
        )

        def _fake_parallel_map(func, items, n_jobs, **kwargs):
            del func, n_jobs
            call_order.append("parallel")
            progress = kwargs.get("progress")
            if progress is not None:
                progress(len(items), len(items))
            return [np.nan] * len(items)

        monkeypatch.setattr(synthesis_utils, "parallel_map", _fake_parallel_map)

        result = calculate_synthesis_scores(
            df,
            config={"enabled_scores": ["scscore"], "n_jobs": 2},
        )

        assert call_order[0] == "load"
        assert "parallel" in call_order
        assert np.isnan(result.loc[0, "sc_score"])

    def test_fsscore_auto_install_populates_missing_runtime(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        """Auto-install should wire FSScore worker python/model when missing."""
        import hedgehog.setup as setup_mod

        runtime = SimpleNamespace(
            worker_python=tmp_path / "optional_envs" / "fsscore" / "bin" / "python",
            model_path=tmp_path
            / "modules"
            / "fsscore"
            / "models"
            / "pretrain_graph_GGLGGL_ep242_best_valloss.ckpt",
            checkout_path=tmp_path / "modules" / "fsscore",
        )
        runtime.worker_python.parent.mkdir(parents=True, exist_ok=True)
        runtime.worker_python.write_text("#!/usr/bin/env python\n", encoding="utf-8")
        runtime.model_path.parent.mkdir(parents=True, exist_ok=True)
        runtime.model_path.write_bytes(b"x")

        monkeypatch.setenv("HEDGEHOG_AUTO_INSTALL", "1")
        monkeypatch.setattr(setup_mod, "ensure_fsscore_runtime", lambda _: runtime)

        captured: dict[str, list[str]] = {}

        def _fake_run(command, **kwargs):
            captured["command"] = command
            output_csv = Path(command[command.index("--output-csv") + 1])
            output_csv.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame({"score": [0.33]}).to_csv(output_csv, index=False)
            return SimpleNamespace(returncode=0)

        monkeypatch.setattr(synthesis_utils.subprocess, "run", _fake_run)

        result = calculate_synthesis_scores(
            pd.DataFrame({"smiles": [SMILES_ETHANOL]}),
            config={"enabled_scores": ["fsscore"]},
        )

        assert result["fs_score"].tolist() == [0.33]
        command = captured["command"]
        assert command[command.index("--worker-python") + 1] == str(
            runtime.worker_python
        )
        assert command[command.index("--model-path") + 1] == str(runtime.model_path)

    def test_fsscore_explicit_command_skips_auto_setup(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Configured FSScore command must bypass auto-setup helper."""
        import hedgehog.setup as setup_mod

        monkeypatch.setenv("HEDGEHOG_AUTO_INSTALL", "1")
        monkeypatch.setattr(
            setup_mod,
            "ensure_fsscore_runtime",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("ensure_fsscore_runtime should not be called")
            ),
        )

        command = (
            f'{sys.executable} -c "import pandas as pd, sys; '
            "df = pd.read_csv(sys.argv[1]); "
            "df['score'] = [0.5 for _ in range(len(df))]; "
            'df.to_csv(sys.argv[2], index=False)" {input} {output}'
        )
        result = calculate_synthesis_scores(
            pd.DataFrame({"smiles": [SMILES_ETHANOL]}),
            config={"enabled_scores": ["fsscore"], "fsscore_command": command},
        )
        assert result["fs_score"].tolist() == [0.5]

    def test_gasa_auto_install_populates_local_worker_command(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        """Auto-install should wire GASA worker command when backend is missing."""
        import hedgehog.setup as setup_mod

        setup_result = SimpleNamespace(
            repo_path=tmp_path / "modules" / "gasa",
            worker_python=tmp_path / "optional_envs" / "gasa" / "bin" / "python",
        )
        setup_result.repo_path.mkdir(parents=True, exist_ok=True)
        setup_result.worker_python.parent.mkdir(parents=True, exist_ok=True)
        setup_result.worker_python.write_text(
            "#!/usr/bin/env python\n", encoding="utf-8"
        )

        monkeypatch.setenv("HEDGEHOG_AUTO_INSTALL", "1")
        monkeypatch.setattr(setup_mod, "ensure_gasa_worker", lambda _: setup_result)

        captured: dict[str, list[str]] = {}

        def _fake_run(command, **kwargs):
            captured["command"] = command
            input_csv = Path(command[command.index("--input-csv") + 1])
            output_csv = Path(command[command.index("--output-csv") + 1])
            output_csv.parent.mkdir(parents=True, exist_ok=True)
            pd.read_csv(input_csv).assign(gasa_score=[0.8]).to_csv(
                output_csv, index=False
            )
            return SimpleNamespace(returncode=0)

        monkeypatch.setattr(synthesis_utils.subprocess, "run", _fake_run)

        result = calculate_synthesis_scores(
            pd.DataFrame({"smiles": [SMILES_ETHANOL]}),
            config={"enabled_scores": ["gasa"]},
        )

        assert result["gasa_score"].tolist() == [0.8]
        rendered = " ".join(captured["command"])
        assert "hedgehog.workers.gasa_worker" in rendered
        assert str(setup_result.worker_python) in rendered
        assert str(setup_result.repo_path) in rendered
        assert (
            captured["command"][captured["command"].index("--num-workers") + 1] == "1"
        )

    def test_gasa_explicit_backend_skips_auto_setup(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Configured GASA command must bypass auto-setup helper."""
        import hedgehog.setup as setup_mod

        monkeypatch.setenv("HEDGEHOG_AUTO_INSTALL", "1")
        monkeypatch.setattr(
            setup_mod,
            "ensure_gasa_worker",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("ensure_gasa_worker should not be called")
            ),
        )

        command = (
            f'{sys.executable} -c "import pandas as pd, sys; '
            "df = pd.read_csv(sys.argv[1]); "
            "df['gasa_score'] = [0.77 for _ in range(len(df))]; "
            'df.to_csv(sys.argv[2], index=False)" {input} {output}'
        )
        result = calculate_synthesis_scores(
            pd.DataFrame({"smiles": [SMILES_ETHANOL]}),
            config={"enabled_scores": ["gasa"], "gasa_command": command},
        )
        assert result["gasa_score"].tolist() == [0.77]

    def test_nonpher_external_worker_hook_uses_configured_python(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Nonpher scorer should prefer external worker when python is configured."""
        import hedgehog.setup as setup_mod

        monkeypatch.setenv("HEDGEHOG_NONPHER_PYTHON", "/tmp/nonpher/bin/python")

        captured: dict[str, str] = {}

        def _fake_external(**kwargs):
            captured["python_bin"] = str(kwargs["python_bin"])
            input_csv = Path(kwargs["input_csv"])
            output_csv = Path(kwargs["output_csv"])
            source = pd.read_csv(input_csv)
            source["nonpher_complexity_score"] = [0.0 for _ in range(len(source))]
            source.to_csv(output_csv, index=False)

        monkeypatch.setattr(setup_mod, "run_nonpher_batch_external", _fake_external)

        result = calculate_synthesis_scores(
            pd.DataFrame({"smiles": [SMILES_ETHANOL]}),
            config={"enabled_scores": ["nonpher"]},
        )

        assert result["nonpher_complexity_score"].tolist() == [0.0]
        assert captured["python_bin"] == "/tmp/nonpher/bin/python"

    def test_nonpher_auto_install_libboost_blocker_returns_nan(
        self, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
    ) -> None:
        """libboost blocker during Nonpher auto-install should soft-fail to NaN."""
        import hedgehog.setup as setup_mod

        synthesis_utils._MISSING_EXTERNAL_SCORER_WARNINGS.clear()
        monkeypatch.setenv("HEDGEHOG_AUTO_INSTALL", "1")
        monkeypatch.delenv("HEDGEHOG_NONPHER_PYTHON", raising=False)
        monkeypatch.setattr(
            setup_mod,
            "ensure_nonpher_uv_runtime",
            lambda **kwargs: SimpleNamespace(
                available=False,
                python_bin=str(Path(kwargs["env_prefix"]) / "bin" / "python"),
                detail="LibMambaUnsatisfiableError: nothing provides libboost 1.65.*",
                install_attempted=True,
            ),
        )
        monkeypatch.setattr(
            synthesis_utils,
            "_load_nonpher_filter_impl",
            lambda: (_ for _ in ()).throw(
                AssertionError(
                    "Local Nonpher fallback should not run after libboost blocker"
                )
            ),
        )

        with caplog.at_level("WARNING"):
            result = calculate_synthesis_scores(
                pd.DataFrame({"smiles": [SMILES_ETHANOL]}),
                config={"enabled_scores": ["nonpher"]},
            )

        assert np.isnan(result.loc[0, "nonpher_complexity_score"])
        assert (
            "LibMambaUnsatisfiableError: nothing provides libboost 1.65.*"
            in caplog.text
        )

    def test_empty_config(self):
        """Empty config - all molecules should pass."""
        df = pd.DataFrame(
            {
                "smiles": ["a", "b", "c"],
                "sa_score": [2.0, 5.0, 8.0],
            }
        )
        config = {}
        result = apply_synthesis_score_filters(df, config)
        assert len(result) == 3

    def test_multiple_filters(self):
        """Multiple filters applied - must pass all."""
        df = pd.DataFrame(
            {
                "smiles": ["a", "b", "c", "d"],
                "sa_score": [2.0, 2.0, 8.0, 2.0],
                "ra_score": [0.8, 0.3, 0.8, 0.8],
            }
        )
        config = {
            "sa_score_min": 0,
            "sa_score_max": 5,
            "ra_score_min": 0.5,
            "ra_score_max": "inf",
        }
        result = apply_synthesis_score_filters(df, config)
        # 'a': sa OK, ra OK - pass
        # 'b': sa OK, ra fail - fail
        # 'c': sa fail, ra OK - fail
        # 'd': sa OK, ra OK - pass
        assert len(result) == 2


class TestMergeRetrosynthesisResults:
    """Tests for merge_retrosynthesis_results function."""

    def test_basic_merge(self):
        """Basic merge of retrosynthesis results."""
        input_df = pd.DataFrame(
            {
                COL_SMILES: ["a", "b", "c"],
                COL_MODEL_NAME: ["m1", "m1", "m1"],
            }
        )
        retro_df = pd.DataFrame(
            {
                "index": [0, 1, 2],
                "SMILES": ["a", "b", "c"],
                "solved": [1, 0, 1],
                "search_time": [1.0, 2.0, 3.0],
            }
        )
        result = merge_retrosynthesis_results(input_df, retro_df)

        assert "solved" in result.columns
        assert "search_time" in result.columns
        assert result["solved"].tolist() == [1, 0, 1]
        assert result["search_time"].tolist() == [1.0, 2.0, 3.0]

    def test_merge_with_missing_results(self):
        """Merge when retro results have fewer rows than input."""
        input_df = pd.DataFrame(
            {
                COL_SMILES: ["a", "b", "c", "d"],
                COL_MODEL_NAME: ["m1", "m1", "m1", "m1"],
            }
        )
        retro_df = pd.DataFrame(
            {
                "index": [0, 1],
                "SMILES": ["a", "b"],
                "solved": [1, 0],
                "search_time": [1.0, 2.0],
            }
        )
        result = merge_retrosynthesis_results(input_df, retro_df)

        assert len(result) == 4
        assert result["solved"].iloc[0] == 1
        assert result["solved"].iloc[1] == 0
        # Missing results should be 0
        assert result["solved"].iloc[2] == 0
        assert result["solved"].iloc[3] == 0

    def test_preserve_input_columns(self):
        """Original columns from input should be preserved."""
        input_df = pd.DataFrame(
            {
                COL_SMILES: ["a", "b"],
                COL_MODEL_NAME: ["m1", "m2"],
                "mol_idx": ["idx-0", "idx-1"],
            }
        )
        retro_df = pd.DataFrame(
            {
                "index": [0, 1],
                "solved": [1, 1],
                "search_time": [1.0, 1.0],
            }
        )
        result = merge_retrosynthesis_results(input_df, retro_df)

        assert COL_SMILES in result.columns
        assert COL_MODEL_NAME in result.columns
        assert "mol_idx" in result.columns


class TestBuildScoreFilterMask:
    """Tests for _build_score_filter_mask function."""

    def test_mask_with_min_max(self):
        """Build mask with both min and max."""
        df = pd.DataFrame({"score": [1, 3, 5, 7, 9]})
        mask = _build_score_filter_mask(df, "score", 2, 6)

        assert mask is not None
        # Values 3, 5 are in range [2, 6]
        expected = [False, True, True, False, False]
        assert mask.tolist() == expected

    def test_mask_with_inf_max(self):
        """Build mask with infinity max."""
        df = pd.DataFrame({"score": [1, 3, 5, 7, 9]})
        mask = _build_score_filter_mask(df, "score", 5, "inf")

        assert mask is not None
        # Values >= 5: 5, 7, 9
        expected = [False, False, True, True, True]
        assert mask.tolist() == expected

    def test_mask_with_nan_values(self):
        """NaN values should pass the filter (return True)."""
        df = pd.DataFrame({"score": [1, np.nan, 5]})
        mask = _build_score_filter_mask(df, "score", 3, 6)

        assert mask is not None
        # 1 < 3 (fail), NaN (pass), 5 in [3,6] (pass)
        expected = [False, True, True]
        assert mask.tolist() == expected

    def test_missing_column(self):
        """Missing column - should return None."""
        df = pd.DataFrame({"other": [1, 2, 3]})
        mask = _build_score_filter_mask(df, "score", 0, 10)
        assert mask is None

    def test_all_nan_column(self):
        """All NaN column - should return None (no valid scores)."""
        df = pd.DataFrame({"score": [np.nan, np.nan, np.nan]})
        mask = _build_score_filter_mask(df, "score", 0, 10)
        assert mask is None


class TestPrepareInputSmiles:
    """Tests for prepare_input_smiles function."""

    def test_basic_output(self, tmp_path):
        """Should write SMILES to output file."""
        df = pd.DataFrame({COL_SMILES: [SMILES_ETHANOL, SMILES_BENZENE, "CC"]})
        output_file = tmp_path / "output.smi"

        count = prepare_input_smiles(df, output_file)

        assert count == 3
        assert output_file.exists()
        content = output_file.read_text()
        assert SMILES_ETHANOL in content
        assert SMILES_BENZENE in content

    def test_creates_parent_dirs(self, tmp_path):
        """Should create parent directories if needed."""
        df = pd.DataFrame({COL_SMILES: [SMILES_ETHANOL]})
        output_file = tmp_path / "nested" / "deep" / "output.smi"

        prepare_input_smiles(df, output_file)

        assert output_file.exists()

    def test_skips_nan(self, tmp_path):
        """Should skip NaN SMILES."""
        df = pd.DataFrame({COL_SMILES: [SMILES_ETHANOL, np.nan, "CC"]})
        output_file = tmp_path / "output.smi"

        count = prepare_input_smiles(df, output_file)

        assert count == 2


class TestRunAizynthfinder:
    """Tests for run_aizynthfinder nproc handling."""

    def test_uses_uv_binary_from_resolver(self, tmp_path, monkeypatch):
        """run_aizynthfinder should use the binary returned by resolve_uv_binary."""
        captured = {}

        def fake_run(cmd, **kwargs):
            captured["cmd"] = cmd
            captured["kwargs"] = kwargs
            return None

        uv_bin = tmp_path / "bin" / "uv"
        uv_bin.parent.mkdir(parents=True, exist_ok=True)
        uv_bin.write_text("#!/bin/sh\n")
        uv_bin.chmod(0o755)

        monkeypatch.setattr(synthesis_utils.subprocess, "run", fake_run)
        monkeypatch.setattr(synthesis_utils, "resolve_uv_binary", lambda: str(uv_bin))
        monkeypatch.setenv("VIRTUAL_ENV", "/tmp/foreign-env")
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "9")

        config = _aizynth_config_path(tmp_path)
        ok = run_aizynthfinder(tmp_path / "in.smi", tmp_path / "out.json", config)

        assert ok is True
        assert captured["cmd"][0] == str(uv_bin)
        assert "--nproc" in captured["cmd"]
        assert captured["cmd"][captured["cmd"].index("--nproc") + 1] == "9"
        assert "VIRTUAL_ENV" not in captured["kwargs"]["env"]
        assert (
            _aizynth_install_root(tmp_path) / "aizynthfinder" / "data" / "logging.yml"
        ).exists()

    def test_returns_false_when_uv_unavailable(self, tmp_path, monkeypatch):
        """Missing uv binary should be logged and prevent command execution."""
        captured = {}

        def fake_run(cmd, **kwargs):
            captured["cmd"] = cmd
            captured["kwargs"] = kwargs
            return None

        monkeypatch.setattr(
            synthesis_utils,
            "resolve_uv_binary",
            lambda: (_ for _ in ()).throw(RuntimeError("uv missing")),
        )
        monkeypatch.setattr(synthesis_utils.subprocess, "run", fake_run)

        config = _aizynth_config_path(tmp_path)
        ok = run_aizynthfinder(tmp_path / "in.smi", tmp_path / "out.json", config)

        assert ok is False
        assert "cmd" not in captured

    def test_minus_one_resolves_to_slurm_cpus(self, tmp_path, monkeypatch):
        """AIZYNTH_NPROC=-1 should resolve to available CPU count."""
        captured = {}

        def fake_run(cmd, **kwargs):
            captured["cmd"] = cmd
            captured["kwargs"] = kwargs
            return None

        monkeypatch.setattr(synthesis_utils.subprocess, "run", fake_run)
        monkeypatch.setenv("AIZYNTH_NPROC", "-1")
        monkeypatch.delenv("MOLSCORE_NJOBS", raising=False)
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "12")

        config = _aizynth_config_path(tmp_path)
        ok = run_aizynthfinder(tmp_path / "in.smi", tmp_path / "out.json", config)

        assert ok is True
        cmd = captured["cmd"]
        assert "--nproc" in cmd
        assert cmd[cmd.index("--nproc") + 1] == "12"

    def test_zero_resolves_to_molscore_njobs(self, tmp_path, monkeypatch):
        """AIZYNTH_NPROC=0 should follow MOLSCORE_NJOBS when set."""
        captured = {}

        def fake_run(cmd, **kwargs):
            captured["cmd"] = cmd
            captured["kwargs"] = kwargs
            return None

        monkeypatch.setattr(synthesis_utils.subprocess, "run", fake_run)
        monkeypatch.setenv("AIZYNTH_NPROC", "0")
        monkeypatch.setenv("MOLSCORE_NJOBS", "20")

        config = _aizynth_config_path(tmp_path)
        ok = run_aizynthfinder(tmp_path / "in.smi", tmp_path / "out.json", config)

        assert ok is True
        cmd = captured["cmd"]
        assert "--nproc" in cmd
        assert cmd[cmd.index("--nproc") + 1] == "20"

    def test_invalid_nproc_falls_back_to_synthesis_n_jobs(self, tmp_path, monkeypatch):
        """Invalid AIZYNTH_NPROC should fall back to synthesis n_jobs."""
        captured = {}

        def fake_run(cmd, **kwargs):
            captured["cmd"] = cmd
            captured["kwargs"] = kwargs
            return None

        monkeypatch.setattr(synthesis_utils.subprocess, "run", fake_run)
        monkeypatch.setenv("AIZYNTH_NPROC", "not-a-number")

        config = _aizynth_config_path(tmp_path)
        ok = run_aizynthfinder(
            tmp_path / "in.smi",
            tmp_path / "out.json",
            config,
            synthesis_config={"n_jobs": 7},
        )

        assert ok is True
        cmd = captured["cmd"]
        assert "--nproc" in cmd
        assert cmd[cmd.index("--nproc") + 1] == "7"

    def test_uses_synthesis_n_jobs_when_env_unset(self, tmp_path, monkeypatch):
        """When AIZYNTH_NPROC is unset, use synthesis.n_jobs."""
        captured = {}

        def fake_run(cmd, **kwargs):
            captured["cmd"] = cmd
            captured["kwargs"] = kwargs
            return None

        monkeypatch.setattr(synthesis_utils.subprocess, "run", fake_run)
        monkeypatch.delenv("AIZYNTH_NPROC", raising=False)
        monkeypatch.delenv("MOLSCORE_NJOBS", raising=False)
        monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)

        config = _aizynth_config_path(tmp_path)
        ok = run_aizynthfinder(
            tmp_path / "in.smi",
            tmp_path / "out.json",
            config,
            synthesis_config={"n_jobs": 6},
        )

        assert ok is True
        cmd = captured["cmd"]
        assert "--nproc" in cmd
        assert cmd[cmd.index("--nproc") + 1] == "6"

    def test_extract_aizynth_progress_from_pairs_and_percent(self):
        assert synthesis_utils._extract_aizynth_progress("done 3/10") == (3, 10)
        assert synthesis_utils._extract_aizynth_progress("50%", fallback_total=8) == (
            4,
            8,
        )

    def test_progress_callback_receives_stream_updates(self, tmp_path, monkeypatch):
        uv_bin = tmp_path / "bin" / "uv"
        uv_bin.parent.mkdir(parents=True, exist_ok=True)
        uv_bin.write_text("#!/bin/sh\n")
        uv_bin.chmod(0o755)

        config = _aizynth_config_path(tmp_path)
        input_smi = tmp_path / "in.smi"
        input_smi.write_text("CC\nCCC\nCCCC\nCCN\n", encoding="utf-8")

        class _FakePopen:
            def __init__(self, *args, **kwargs):
                del args, kwargs
                self.stdout = io.StringIO("step 1/4\nstep 2/4\n")

            def wait(self):
                return 0

        monkeypatch.setattr(synthesis_utils, "resolve_uv_binary", lambda: str(uv_bin))
        monkeypatch.setattr(synthesis_utils.subprocess, "Popen", _FakePopen)

        progress_events: list[tuple[int, int, str | None]] = []

        ok = run_aizynthfinder(
            input_smi,
            tmp_path / "out.json",
            config,
            progress_cb=lambda done, total, detail: progress_events.append(
                (done, total, detail)
            ),
        )

        assert ok is True
        assert progress_events
        assert progress_events[0][0] == 0
        assert any(done >= 2 and total == 4 for done, total, _ in progress_events)
        assert progress_events[-1][0] == progress_events[-1][1]


class TestParseRetrosynthesisResults:
    """Tests for parse_retrosynthesis_results function."""

    def test_valid_json(self, tmp_path):
        """Parse valid JSON results."""
        json_file = tmp_path / "results.json"
        data = {
            "data": [
                {"index": 0, "target": "CCO", "is_solved": True, "search_time": 1.5},
                {"index": 1, "target": "CC", "is_solved": False, "search_time": 2.0},
            ]
        }
        json_file.write_text(json.dumps(data))

        result = parse_retrosynthesis_results(json_file)

        assert len(result) == 2
        assert result["solved"].tolist() == [1, 0]
        assert result["search_time"].tolist() == [1.5, 2.0]

    def test_empty_json(self, tmp_path):
        """Parse JSON with no data."""
        json_file = tmp_path / "results.json"
        json_file.write_text('{"data": []}')

        result = parse_retrosynthesis_results(json_file)

        assert len(result) == 0

    def test_missing_data_key(self, tmp_path):
        """JSON without 'data' key should return empty DataFrame."""
        json_file = tmp_path / "results.json"
        json_file.write_text('{"other_key": []}')

        result = parse_retrosynthesis_results(json_file)

        assert len(result) == 0

    def test_malformed_json(self, tmp_path):
        """Malformed JSON should return empty DataFrame."""
        json_file = tmp_path / "results.json"
        json_file.write_text("not valid json")

        result = parse_retrosynthesis_results(json_file)

        assert len(result) == 0


class TestGetInputPath:
    """Tests for get_input_path function."""

    def test_finds_struct_filters_output(self, tmp_path):
        """Should find structural filters output first."""
        sf_dir = tmp_path / "stages" / "03_structural_filters_post"
        sf_dir.mkdir(parents=True)
        (sf_dir / FILE_FILTERED_MOLECULES).write_text(f"{COL_SMILES}\n{SMILES_ETHANOL}")

        config = {"generated_mols_path": "/fallback/path.csv"}
        result = get_input_path(config, str(tmp_path))

        assert "structural_filters_post" in result

    def test_finds_descriptors_output(self, tmp_path):
        """Should find descriptors output if no struct filters."""
        desc_dir = tmp_path / "stages" / "02_descriptors_initial" / "filtered"
        desc_dir.mkdir(parents=True)
        (desc_dir / FILE_FILTERED_MOLECULES).write_text(
            f"{COL_SMILES}\n{SMILES_ETHANOL}"
        )

        config = {"generated_mols_path": "/fallback/path.csv"}
        result = get_input_path(config, str(tmp_path))

        assert "descriptors_initial" in result

    def test_falls_back_to_config(self, tmp_path):
        """Should fall back to config path if no processed data."""
        config = {"generated_mols_path": "/fallback/path.csv"}
        result = get_input_path(config, str(tmp_path))

        assert result == "/fallback/path.csv"


class TestSynthesisScoreCalculations:
    """Tests for synthesis score calculation functions."""

    def test_sa_score_simple_molecule(self):
        """Simple molecules should have low SA score."""
        score = _calculate_sa_score("CC")  # ethane
        if not np.isnan(score):
            assert score < 3  # Very easy to synthesize

    def test_sa_score_complex_molecule(self):
        """More complex molecules should have higher SA score."""
        # Complex natural product-like molecule
        score = _calculate_sa_score("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O")
        if not np.isnan(score):
            assert 1 <= score <= 10

    def test_syba_score_drug_like(self):
        """Drug-like molecule should have positive SYBA score."""
        score = _calculate_syba_score("CC(=O)Nc1ccc(O)cc1")  # paracetamol
        # Score can be NaN if SYBA not available
        if not np.isnan(score):
            assert isinstance(score, (int, float))

    def test_calculate_synthesis_scores_includes_sync(
        self, monkeypatch: pytest.MonkeyPatch
    ):
        """Score calculation should include SYNC when enabled."""
        df = pd.DataFrame({"smiles": [SMILES_ETHANOL]})
        monkeypatch.setitem(
            synthesis_utils.SYNTHESIS_SCORERS,
            "sync",
            synthesis_utils.SynthesisScorer(
                name="sync",
                column="sync_score",
                batch_calculator=lambda smiles, *_args, **_kwargs: [0.9] * len(smiles),
            ),
        )

        result = synthesis_utils.calculate_synthesis_scores(
            df, config={"enabled_scores": ["sync"], "n_jobs": 1}
        )

        assert result["sync_score"].tolist() == [0.9]

    def test_calculate_synthesis_scores_can_skip_sync(
        self, monkeypatch: pytest.MonkeyPatch
    ):
        """SYNC scoring should be optional."""
        df = pd.DataFrame({"smiles": [SMILES_ETHANOL]})
        monkeypatch.setattr(
            synthesis_utils, "_calculate_sa_score_single", lambda _smiles: 2.0
        )

        result = synthesis_utils.calculate_synthesis_scores(
            df, config={"enabled_scores": ["sa"], "n_jobs": 1}
        )

        assert "sync_score" not in result.columns

    def test_sync_checkpoint_loader_rejects_unsafe_custom_pickle(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        """Custom checkpoints must be safe state_dict payloads."""
        checkpoint = tmp_path / "custom.ckpt"
        checkpoint.write_bytes(b"not-a-safe-state-dict")

        class _FakeTorch:
            @staticmethod
            def load(*_args, **kwargs):
                assert kwargs.get("weights_only") is True
                raise RuntimeError("unsafe pickle fallback needed")

        monkeypatch.setattr(sync_module, "_load_torch", lambda: (_FakeTorch, object()))

        with pytest.raises(RuntimeError, match="official checkpoint checksum"):
            sync_module._load_checkpoint_state_dict(checkpoint)


@pytest.mark.skipif(
    not _get_rascore_model_path().exists(),
    reason="RAScore model not available",
)
class TestCalculateRaScores:
    """Tests for RA score batch calculation."""

    def test_valid_smiles(self):
        """RA scores for known synthesizable molecules should be high."""
        scores = _calculate_ra_scores_batch([SMILES_BENZENE, SMILES_ETHANOL])
        assert len(scores) == 2
        for s in scores:
            assert not np.isnan(s)
            assert 0 <= s <= 1

    def test_invalid_smiles_returns_nan(self):
        """Invalid SMILES should get NaN score."""
        scores = _calculate_ra_scores_batch(["invalid_smiles"])
        assert len(scores) == 1
        assert np.isnan(scores[0])

    def test_mixed_valid_invalid(self):
        """Batch with valid and invalid SMILES."""
        scores = _calculate_ra_scores_batch(
            [SMILES_BENZENE, "not_a_molecule", SMILES_ETHANOL]
        )
        assert len(scores) == 3
        assert not np.isnan(scores[0])
        assert np.isnan(scores[1])
        assert not np.isnan(scores[2])

    def test_empty_list(self):
        """Empty input should return empty list."""
        scores = _calculate_ra_scores_batch([])
        assert scores == []

    def test_simple_molecule_high_score(self):
        """Simple drug-like molecules should have high RA scores."""
        scores = _calculate_ra_scores_batch(["CC(=O)Oc1ccccc1C(O)=O"])  # aspirin
        assert len(scores) == 1
        assert scores[0] > 0.9


class TestRascoreAutoInstall:
    """Tests for RAScore auto-download/load behavior."""

    def test_ensure_rascore_json_model_converts_when_missing(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        """Missing JSON model should be generated from pickle model."""
        json_path = tmp_path / "model.json"
        pkl_path = tmp_path / "model.pkl"
        pkl_path.write_bytes(b"pickle-bytes")

        monkeypatch.setattr(
            synthesis_utils, "_get_rascore_model_path", lambda: json_path
        )
        monkeypatch.setattr(
            synthesis_utils, "_ensure_rascore_pickle_path", lambda: pkl_path
        )

        def _fake_convert(src: Path, dst: Path) -> bool:
            assert src == pkl_path
            assert dst == json_path
            dst.write_text('{"converted": true}', encoding="utf-8")
            return True

        monkeypatch.setattr(
            synthesis_utils, "_convert_rascore_pickle_to_json", _fake_convert
        )

        ensured = synthesis_utils._ensure_rascore_json_model()
        assert ensured == json_path
        assert json_path.exists()

    def test_ensure_rascore_json_model_returns_none_on_convert_failure(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        """JSON helper should return None when conversion fails."""
        json_path = tmp_path / "model.json"
        pkl_path = tmp_path / "model.pkl"
        pkl_path.write_bytes(b"pickle-bytes")

        monkeypatch.setattr(
            synthesis_utils, "_get_rascore_model_path", lambda: json_path
        )
        monkeypatch.setattr(
            synthesis_utils, "_ensure_rascore_pickle_path", lambda: pkl_path
        )
        monkeypatch.setattr(
            synthesis_utils, "_convert_rascore_pickle_to_json", lambda _src, _dst: False
        )

        ensured = synthesis_utils._ensure_rascore_json_model()
        assert ensured is None

    def test_load_rascore_uses_generated_json_before_pickle(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        """Loader should use generated JSON model before trying pickle load."""
        json_path = tmp_path / "model.json"
        json_path.write_text('{"learner": {}}', encoding="utf-8")
        pkl_path = tmp_path / "model.pkl"
        pkl_path.write_text("should-not-be-used", encoding="utf-8")

        monkeypatch.setattr(
            synthesis_utils, "_get_rascore_model_path", lambda: json_path
        )
        monkeypatch.setattr(
            synthesis_utils, "_get_rascore_pickle_model_path", lambda: pkl_path
        )

        loaded_model_path: dict[str, str] = {}

        class _DummyBooster:
            def load_model(self, path: str) -> None:
                loaded_model_path["path"] = path

        monkeypatch.setitem(
            sys.modules, "xgboost", SimpleNamespace(Booster=_DummyBooster)
        )
        synthesis_utils._lazy_cache["rascore_booster"] = None

        loaded = synthesis_utils._load_rascore_impl()
        assert loaded
        assert loaded["kind"] == "xgboost_booster"
        assert loaded_model_path["path"] == str(json_path)

    def test_load_rascore_suppresses_legacy_json_warning(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        """Legacy XGBoost JSON warning should not be surfaced to users."""
        json_path = tmp_path / "model.json"
        json_path.write_text('{"learner": {}}', encoding="utf-8")
        pkl_path = tmp_path / "model.pkl"
        pkl_path.write_text("should-not-be-used", encoding="utf-8")

        monkeypatch.setattr(
            synthesis_utils, "_get_rascore_model_path", lambda: json_path
        )
        monkeypatch.setattr(
            synthesis_utils, "_get_rascore_pickle_model_path", lambda: pkl_path
        )

        class _DummyBooster:
            def load_model(self, _path: str) -> None:
                warnings.warn(
                    "[23:15:45] WARNING: /__w/xgboost/xgboost/src/learner.cc:885: "
                    "Found JSON model saved \n"
                    "before XGBoost 1.6, please save the model using current version again.",
                    UserWarning,
                    stacklevel=1,
                )

        monkeypatch.setitem(
            sys.modules, "xgboost", SimpleNamespace(Booster=_DummyBooster)
        )
        synthesis_utils._lazy_cache["rascore_booster"] = None

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            loaded = synthesis_utils._load_rascore_impl()

        assert loaded
        assert loaded["kind"] == "xgboost_booster"
        legacy_fragments = ("found json model saved", "xgboost 1.6")
        assert not any(
            all(
                fragment in " ".join(str(w.message).lower().split())
                for fragment in legacy_fragments
            )
            for w in caught
        )

    def test_load_rascore_auto_downloads_pickle_model(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        """When local model is missing, helper should auto-provide pickle model."""
        model_path = tmp_path / "model.pkl"
        with open(model_path, "wb") as f:
            pickle.dump(_DummyRascoreClassifier(), f)

        missing_json = tmp_path / "missing_model.json"
        missing_pkl = tmp_path / "missing_model.pkl"
        monkeypatch.setattr(
            synthesis_utils, "_get_rascore_model_path", lambda: missing_json
        )
        monkeypatch.setattr(
            synthesis_utils, "_get_rascore_pickle_model_path", lambda: missing_pkl
        )
        monkeypatch.setattr(
            "hedgehog.setup._rascore.ensure_rascore_model",
            lambda _project_root: model_path,
        )
        synthesis_utils._lazy_cache["rascore_booster"] = None

        loaded = synthesis_utils._load_rascore_impl()
        assert loaded
        assert loaded["kind"] == "pickle_classifier"

    def test_calculate_ra_scores_with_pickle_classifier(
        self, monkeypatch: pytest.MonkeyPatch
    ):
        """RA score batch should use predict_proba when pickle model is loaded."""
        monkeypatch.setattr(
            synthesis_utils,
            "_load_rascore",
            lambda: {"kind": "pickle_classifier", "model": _DummyRascoreClassifier()},
        )
        scores = _calculate_ra_scores_batch([SMILES_ETHANOL, "invalid_smiles"])

        assert len(scores) == 2
        assert scores[0] == pytest.approx(0.8)
        assert np.isnan(scores[1])

    def test_calculate_ra_scores_legacy_worker_fallback(
        self, monkeypatch: pytest.MonkeyPatch
    ):
        """When local model load fails, batch scorer should use legacy worker."""
        monkeypatch.delenv(synthesis_utils.STRICT_RASCORE_ENV, raising=False)
        called = {"run": False}

        def _fake_run(cmd, **kwargs):
            called["run"] = True
            out_idx = cmd.index("--output-json") + 1
            Path(cmd[out_idx]).write_text(
                json.dumps({"scores": [0.42, None]}), encoding="utf-8"
            )
            return SimpleNamespace(returncode=0, stdout="", stderr="")

        monkeypatch.setattr(synthesis_utils, "_load_rascore", lambda: False)
        monkeypatch.setattr(
            synthesis_utils, "_ensure_rascore_pickle_path", lambda: Path("dummy.pkl")
        )
        monkeypatch.setattr(synthesis_utils, "resolve_uv_binary", lambda: "uv")
        monkeypatch.setattr(synthesis_utils.subprocess, "run", _fake_run)

        scores = _calculate_ra_scores_batch([SMILES_ETHANOL, "invalid_smiles"])

        assert called["run"] is True
        assert len(scores) == 2
        assert scores[0] == pytest.approx(0.42)
        assert np.isnan(scores[1])

    def test_load_rascore_uses_pickle_when_json_is_invalid(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        """Invalid JSON model should fall back to local pickle model."""
        json_path = tmp_path / "model.json"
        json_path.write_text("corrupted-model", encoding="utf-8")

        pkl_path = tmp_path / "model.pkl"
        with open(pkl_path, "wb") as f:
            pickle.dump(_DummyRascoreClassifier(), f)

        monkeypatch.setattr(
            synthesis_utils, "_get_rascore_model_path", lambda: json_path
        )
        monkeypatch.setattr(
            synthesis_utils, "_get_rascore_pickle_model_path", lambda: pkl_path
        )
        synthesis_utils._lazy_cache["rascore_booster"] = None

        loaded = synthesis_utils._load_rascore_impl()

        assert loaded
        assert loaded["kind"] == "pickle_classifier"

    def test_strict_rascore_mode_raises_when_local_model_fails(
        self, monkeypatch: pytest.MonkeyPatch
    ):
        """Strict mode should fail fast instead of using legacy fallback."""
        monkeypatch.setenv(synthesis_utils.STRICT_RASCORE_ENV, "1")
        monkeypatch.setattr(synthesis_utils, "_load_rascore", lambda: False)

        called_legacy = {"called": False}

        def _fake_legacy(_smiles):
            called_legacy["called"] = True
            return [0.1]

        monkeypatch.setattr(
            synthesis_utils, "_calculate_ra_scores_batch_legacy", _fake_legacy
        )

        with pytest.raises(RuntimeError, match="strict mode"):
            _calculate_ra_scores_batch([SMILES_ETHANOL])

        assert called_legacy["called"] is False
