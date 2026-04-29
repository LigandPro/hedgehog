"""Helpers for validating optional Nonpher runtime availability."""

from __future__ import annotations

import os
import shlex
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from rdkit import Chem

from hedgehog.setup._download import resolve_uv_binary

NONPHER_PYTHON_ENV_VAR = "HEDGEHOG_NONPHER_PYTHON"
NONPHER_FALLBACK_PYTHON_ENV_VAR = "HEDGEHOG_NONPHER_FALLBACK_PYTHON"
DEFAULT_NONPHER_GIT_URL = "git+https://github.com/lich-uct/nonpher.git"
DEFAULT_MOLPHER_GIT_URL = "git+https://github.com/lich-uct/molpher-lib.git"
DEFAULT_NONPHER_FALLBACK_PYTHON = (
    "/mnt/ligandpro/shared_storage/data/nikolenko/"
    "hedgehog_optional_envs/nonpher-hybrid-py38-v2/bin/python"
)


@dataclass(frozen=True)
class NonpherCheckResult:
    """Result of checking whether Nonpher runtime is available."""

    available: bool
    detail: str


@dataclass(frozen=True)
class NonpherEnsureResult:
    """Result of ensuring optional Nonpher runtime in isolated env."""

    available: bool
    python_bin: str | None
    detail: str
    install_attempted: bool


def _short_exc_message(exc: Exception) -> str:
    raw = str(exc).strip()
    if not raw:
        return exc.__class__.__name__
    return raw.splitlines()[0]


def _first_nonempty_line(text: str) -> str:
    for line in text.splitlines():
        stripped = line.strip()
        if stripped:
            return stripped
    return ""


def _first_matching_line(text: str, fragments: tuple[str, ...]) -> str:
    lowered_fragments = tuple(fragment.lower() for fragment in fragments)
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        lowered = stripped.lower()
        if any(fragment in lowered for fragment in lowered_fragments):
            return stripped
    return ""


def resolve_nonpher_python(python_bin: str | None = None) -> str | None:
    """Resolve external Nonpher interpreter from arg or env var."""
    if python_bin and python_bin.strip():
        return python_bin.strip()
    env_value = os.getenv(NONPHER_PYTHON_ENV_VAR, "").strip()
    return env_value or None


def resolve_nonpher_fallback_python(python_bin: str | None = None) -> str | None:
    """Resolve a validated fallback Nonpher interpreter, if one is available."""
    candidates = (
        python_bin,
        os.getenv(NONPHER_FALLBACK_PYTHON_ENV_VAR, "").strip(),
        DEFAULT_NONPHER_FALLBACK_PYTHON,
    )
    for candidate in candidates:
        if not candidate:
            continue
        expanded = Path(candidate).expanduser()
        if expanded.exists():
            return str(expanded)
    return None


def _check_nonpher_fallback_python(
    fallback_python: str | None,
    probe_smiles: str,
    detail_prefix: str,
) -> NonpherEnsureResult | None:
    resolved_fallback = resolve_nonpher_fallback_python(fallback_python)
    if not resolved_fallback:
        return None

    checked = _check_nonpher_external_python(resolved_fallback, probe_smiles)
    if not checked.available:
        return None

    detail = checked.detail
    if detail_prefix:
        detail = f"{detail_prefix}; using fallback runtime. {checked.detail}"
    return NonpherEnsureResult(
        available=True,
        python_bin=resolved_fallback,
        detail=detail,
        install_attempted=False,
    )


def create_nonpher_complexity_filter() -> Any:
    """Create a Nonpher complexity filter instance in the current environment."""
    from nonpher import nonpher as nonpher_mod

    return nonpher_mod.ComplexityFilter()


def _probe_complexity_filter(complexity_filter: Any, probe_smiles: str) -> bool:
    """Probe a filter instance on one molecule to verify runtime wiring."""
    mol = Chem.MolFromSmiles(probe_smiles)
    if mol is None:
        raise ValueError(f"Invalid probe SMILES: {probe_smiles}")
    return bool(complexity_filter.isTooComplex(mol))


def _check_nonpher_current_env(probe_smiles: str) -> NonpherCheckResult:
    try:
        complexity_filter = create_nonpher_complexity_filter()
        _probe_complexity_filter(complexity_filter, probe_smiles)
        return NonpherCheckResult(
            available=True,
            detail="Nonpher runtime probe succeeded in the current environment.",
        )
    except Exception as exc:
        return NonpherCheckResult(
            available=False,
            detail=f"{exc.__class__.__name__}: {_short_exc_message(exc)}",
        )


def _check_nonpher_external_python(
    python_bin: str,
    probe_smiles: str,
) -> NonpherCheckResult:
    script = (
        "import sys\n"
        "from rdkit import Chem\n"
        "from nonpher import nonpher as nonpher_mod\n"
        "probe = sys.argv[1]\n"
        "mol = Chem.MolFromSmiles(probe)\n"
        "if mol is None:\n"
        "    raise ValueError(f'Invalid probe SMILES: {probe}')\n"
        "flt = nonpher_mod.ComplexityFilter()\n"
        "_ = bool(flt.isTooComplex(mol))\n"
        "print('ok')\n"
    )

    try:
        completed = subprocess.run(
            [python_bin, "-c", script, probe_smiles],
            shell=False,
            capture_output=True,
            text=True,
            check=False,
        )
    except Exception as exc:
        return NonpherCheckResult(
            available=False,
            detail=(
                f"Unable to run interpreter '{python_bin}': {_short_exc_message(exc)}"
            ),
        )

    if completed.returncode == 0 and completed.stdout.strip().endswith("ok"):
        return NonpherCheckResult(
            available=True,
            detail=(
                "Nonpher runtime probe succeeded via external interpreter "
                f"{python_bin}."
            ),
        )

    detail = completed.stderr.strip() or completed.stdout.strip() or "probe failed"
    return NonpherCheckResult(
        available=False,
        detail=f"External probe failed for '{python_bin}': {detail.splitlines()[0]}",
    )


def check_nonpher_runtime(
    *,
    python_bin: str | None = None,
    probe_smiles: str = "CCO",
) -> NonpherCheckResult:
    """Check whether Nonpher runtime is usable in current or external Python env."""
    resolved_python = resolve_nonpher_python(python_bin)
    if resolved_python:
        return _check_nonpher_external_python(resolved_python, probe_smiles)
    return _check_nonpher_current_env(probe_smiles)


def _nonpher_worker_path() -> Path:
    return Path(__file__).resolve().parent.parent / "workers" / "nonpher_worker.py"


def _run_command(
    command: list[str], *, timeout: int = 3600
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        shell=False,
        capture_output=True,
        text=True,
        check=False,
        timeout=timeout,
    )


def _format_command_failure(
    command: list[str], completed: subprocess.CompletedProcess[str]
) -> str:
    detail = (
        _first_matching_line(
            completed.stderr,
            (
                "libboost",
                "cannot find -lmolpher",
                "no module named 'molpher'",
                "no module named 'pkg_resources'",
            ),
        )
        or _first_nonempty_line(completed.stderr)
        or _first_nonempty_line(completed.stdout)
    )
    if not detail:
        detail = f"exit code {completed.returncode}"
    return f"{shlex.join(command)} failed: {detail}"


def ensure_nonpher_external_runtime(
    *,
    python_bin: str | None = None,
    env_prefix: str | None = None,
    conda_bin: str = "conda",
    probe_smiles: str = "CCO",
    timeout: int = 3600,
    nonpher_git_url: str = DEFAULT_NONPHER_GIT_URL,
) -> NonpherEnsureResult:
    """Ensure isolated Nonpher runtime and return probe status + interpreter path."""
    resolved_python = resolve_nonpher_python(python_bin)
    if resolved_python:
        existing = _check_nonpher_external_python(resolved_python, probe_smiles)
        if existing.available:
            return NonpherEnsureResult(
                available=True,
                python_bin=resolved_python,
                detail=existing.detail,
                install_attempted=False,
            )
        if env_prefix is None:
            return NonpherEnsureResult(
                available=False,
                python_bin=resolved_python,
                detail=existing.detail,
                install_attempted=False,
            )

    if env_prefix is None:
        return NonpherEnsureResult(
            available=False,
            python_bin=resolved_python,
            detail=(
                f"Nonpher interpreter is not configured. Set {NONPHER_PYTHON_ENV_VAR} "
                "or pass python_bin/env_prefix."
            ),
            install_attempted=False,
        )

    prefix = Path(env_prefix).expanduser()
    prefix_python = str(prefix / "bin" / "python")

    create_command = [
        conda_bin,
        "create",
        "-y",
        "--override-channels",
        "--prefix",
        str(prefix),
        "-c",
        "lich",
        "-c",
        "conda-forge",
        "python=3.8",
        "molpher-lib=0.0.0b3",
        "rdkit=2020.09.1",
        "boost-cpp=1.74.0",
    ]
    created = _run_command(create_command, timeout=timeout)
    if created.returncode != 0:
        return NonpherEnsureResult(
            available=False,
            python_bin=prefix_python,
            detail=_format_command_failure(create_command, created),
            install_attempted=True,
        )

    install_nonpher_cmd = [
        prefix_python,
        "-m",
        "pip",
        "install",
        "--no-deps",
        nonpher_git_url,
    ]
    installed = _run_command(install_nonpher_cmd, timeout=timeout)
    if installed.returncode != 0:
        return NonpherEnsureResult(
            available=False,
            python_bin=prefix_python,
            detail=_format_command_failure(install_nonpher_cmd, installed),
            install_attempted=True,
        )

    checked = _check_nonpher_external_python(prefix_python, probe_smiles)
    return NonpherEnsureResult(
        available=checked.available,
        python_bin=prefix_python,
        detail=checked.detail,
        install_attempted=True,
    )


def ensure_nonpher_uv_runtime(
    *,
    env_prefix: str,
    probe_smiles: str = "CCO",
    timeout: int = 3600,
    nonpher_git_url: str = DEFAULT_NONPHER_GIT_URL,
    molpher_git_url: str = DEFAULT_MOLPHER_GIT_URL,
    fallback_python: str | None = None,
) -> NonpherEnsureResult:
    """Attempt isolated uv-only Nonpher bootstrap in ``env_prefix``."""
    prefix = Path(env_prefix).expanduser()
    prefix_python = str(prefix / "bin" / "python")
    uv_marker = prefix / "pyvenv.cfg"

    if prefix.exists() and not uv_marker.exists():
        return NonpherEnsureResult(
            available=False,
            python_bin=prefix_python,
            detail=(
                f"Refusing to overwrite non-virtualenv path: {prefix}. "
                "Choose an empty env_prefix or remove the path explicitly."
            ),
            install_attempted=False,
        )

    if Path(prefix_python).exists() and uv_marker.exists():
        existing = _check_nonpher_external_python(prefix_python, probe_smiles)
        if existing.available:
            return NonpherEnsureResult(
                available=True,
                python_bin=prefix_python,
                detail=existing.detail,
                install_attempted=False,
            )
        fallback = _check_nonpher_fallback_python(
            fallback_python,
            probe_smiles,
            f"Existing uv runtime is unavailable ({existing.detail})",
        )
        if fallback is not None:
            return fallback

    try:
        uv_bin = resolve_uv_binary()
    except Exception as exc:
        return NonpherEnsureResult(
            available=False,
            python_bin=prefix_python,
            detail=f"Unable to resolve uv binary: {_short_exc_message(exc)}",
            install_attempted=False,
        )

    commands: list[list[str]] = [
        [uv_bin, "venv", "--python", "3.10", str(prefix)],
        [uv_bin, "pip", "install", "--python", prefix_python, "numpy<2"],
        [uv_bin, "pip", "install", "--python", prefix_python, "rdkit-pypi==2022.9.5"],
        [uv_bin, "pip", "install", "--python", prefix_python, "setuptools<81", "wheel"],
        [
            uv_bin,
            "pip",
            "install",
            "--python",
            prefix_python,
            "--no-deps",
            nonpher_git_url,
        ],
        [
            uv_bin,
            "pip",
            "install",
            "--python",
            prefix_python,
            "--no-build-isolation",
            molpher_git_url,
        ],
    ]

    for command in commands:
        completed = _run_command(command, timeout=timeout)
        if completed.returncode != 0:
            failure_detail = _format_command_failure(command, completed)
            fallback = _check_nonpher_fallback_python(
                fallback_python,
                probe_smiles,
                f"uv-only bootstrap failed ({failure_detail})",
            )
            if fallback is not None:
                return NonpherEnsureResult(
                    available=fallback.available,
                    python_bin=fallback.python_bin,
                    detail=fallback.detail,
                    install_attempted=True,
                )
            return NonpherEnsureResult(
                available=False,
                python_bin=prefix_python,
                detail=failure_detail,
                install_attempted=True,
            )

    checked = _check_nonpher_external_python(prefix_python, probe_smiles)
    if not checked.available:
        fallback = _check_nonpher_fallback_python(
            fallback_python,
            probe_smiles,
            f"uv-only runtime probe failed ({checked.detail})",
        )
        if fallback is not None:
            return NonpherEnsureResult(
                available=fallback.available,
                python_bin=fallback.python_bin,
                detail=fallback.detail,
                install_attempted=True,
            )
    return NonpherEnsureResult(
        available=checked.available,
        python_bin=prefix_python,
        detail=checked.detail,
        install_attempted=True,
    )


def nonpher_batch_worker_command(
    *,
    worker_python: str,
    input_csv: str,
    output_csv: str,
    smiles_column: str = "smiles",
    score_column: str = "nonpher_complexity_score",
) -> list[str]:
    """Build external worker command for batch Nonpher scoring."""
    return [
        worker_python,
        str(_nonpher_worker_path()),
        "--input-csv",
        input_csv,
        "--output-csv",
        output_csv,
        "--smiles-column",
        smiles_column,
        "--score-column",
        score_column,
    ]


def run_nonpher_batch_external(
    *,
    input_csv: str,
    output_csv: str,
    smiles_column: str = "smiles",
    score_column: str = "nonpher_complexity_score",
    python_bin: str | None = None,
    timeout: int = 1800,
) -> None:
    """Run batch Nonpher scoring in configured external interpreter."""
    resolved_python = resolve_nonpher_python(python_bin)
    if not resolved_python:
        raise RuntimeError(
            f"External Nonpher Python is not configured. Set {NONPHER_PYTHON_ENV_VAR} "
            "or pass python_bin."
        )

    command = nonpher_batch_worker_command(
        worker_python=resolved_python,
        input_csv=input_csv,
        output_csv=output_csv,
        smiles_column=smiles_column,
        score_column=score_column,
    )
    completed = _run_command(command, timeout=timeout)
    if completed.returncode != 0:
        detail = _first_nonempty_line(completed.stderr) or _first_nonempty_line(
            completed.stdout
        )
        if not detail:
            detail = f"exit code {completed.returncode}"
        raise RuntimeError(
            f"Nonpher external batch scoring failed for '{resolved_python}': {detail}"
        )


def nonpher_lobachevsky_setup_commands(
    *,
    optional_env_root: str = "~/work/hedgehog_optional_envs",
    project_dir: str = "~/work/Projects/hedgehog",
    fallback_python: str = (
        "/mnt/ligandpro/shared_storage/data/nikolenko/"
        "hedgehog_optional_envs/nonpher-hybrid-py38-v2/bin/python"
    ),
) -> list[str]:
    """Return uv-first Nonpher setup/validation steps for Linux hosts."""
    uv_python = "$HEDGEHOG_OPTIONAL_ENV_ROOT/nonpher/bin/python"
    return [
        "ssh lobachevsky",
        f"cd {project_dir}",
        f"export HEDGEHOG_OPTIONAL_ENV_ROOT={optional_env_root}",
        'mkdir -p "$HEDGEHOG_OPTIONAL_ENV_ROOT"',
        (
            "uv run python - <<'PY'\n"
            "from hedgehog.setup import ensure_nonpher_uv_runtime\n"
            "result = ensure_nonpher_uv_runtime(env_prefix='$HEDGEHOG_OPTIONAL_ENV_ROOT/nonpher')\n"
            "print(result)\n"
            "PY"
        ),
        f"uv run hedgehog setup nonpher-check --python {uv_python}",
        "# If uv-only bootstrap is blocked (for example linker/libboost), use a validated external runtime:",
        f"export {NONPHER_PYTHON_ENV_VAR}={fallback_python}",
        "uv run hedgehog setup nonpher-check",
        f"uv run hedgehog setup nonpher-check --python {fallback_python}",
    ]
