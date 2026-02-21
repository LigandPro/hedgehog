import os
import shutil
from pathlib import Path

from hedgehog.configs.logger import logger


def _validate_optional_tool_path(tool_path, tool_label):
    """Validate optional external tool path and return usable value or None."""
    if not tool_path:
        return None

    path = Path(str(tool_path))
    if path.exists():
        if not path.is_file() or not os.access(path, os.X_OK):
            logger.warning(
                "%s is not executable: %s. Falling back to built-in behavior.",
                tool_label,
                tool_path,
            )
            return None
        return str(path)

    resolved = shutil.which(str(tool_path))
    if resolved:
        return resolved

    logger.warning(
        "%s not found: %s. Falling back to built-in behavior.",
        tool_label,
        tool_path,
    )
    return None


def _is_real_binary(path: str) -> bool:
    """Check if a file is a real compiled binary (ELF), not a script wrapper."""
    try:
        with open(path, "rb") as f:
            header = f.read(4)
        return header == b"\x7fELF"
    except OSError:
        return False


def _resolve_docking_binary(config_path: str, tool_name: str) -> str:
    """Resolve a docking binary path from config or PATH.

    Args:
        config_path: Path from config (absolute path or bare tool name).
        tool_name: Tool name for PATH lookup (e.g. 'smina', 'gnina').

    Returns:
        Resolved absolute path to the binary.

    Raises:
        FileNotFoundError: If the binary cannot be found.
    """
    if os.path.isabs(config_path) and os.path.isfile(config_path):
        return config_path

    found = shutil.which(tool_name)
    if found and _is_real_binary(found):
        return found
    elif found:
        logger.debug(
            "%s found at %s but is a script wrapper, not a real binary — skipping",
            tool_name,
            found,
        )

    if tool_name == "gnina":
        from hedgehog.setup import ensure_gnina

        try:
            return ensure_gnina()
        except RuntimeError as exc:
            raise FileNotFoundError(str(exc)) from exc

    raise FileNotFoundError(
        f"Docking binary '{tool_name}' not found. "
        f"Provide absolute path in config or ensure it's on PATH."
    )


def _resolve_path(path, base_dir):
    """Resolve path to absolute, using base_dir if relative."""
    path_obj = Path(path)
    if path_obj.is_absolute():
        return str(path_obj.resolve())
    return str((base_dir / path).resolve())


def _resolve_autobox_path(autobox_ligand, project_root):
    """Resolve autobox_ligand path to absolute."""
    path = Path(autobox_ligand)
    if path.is_absolute():
        return path if path.exists() else None

    candidate = (project_root / autobox_ligand).resolve()
    if candidate.exists():
        return candidate

    if "data/" in autobox_ligand:
        data_path = autobox_ligand[autobox_ligand.find("data/") :]
        candidate = (project_root / data_path).resolve()
        if candidate.exists():
            return candidate

    return None


def _get_gnina_environment(cfg, base_folder):
    """Get GNINA activation command and LD_LIBRARY_PATH."""
    gnina_config = cfg.get("gnina_config", {})
    env_path = gnina_config.get("env_path") or cfg.get("gnina_env_path")

    gnina_activate = cfg.get("gnina_activate")
    if not gnina_activate and env_path:
        conda_sh = cfg.get("conda_sh")
        if not conda_sh:
            # Auto-detect conda.sh from common installation paths
            for candidate in ("miniforge", "miniconda3", "mambaforge", "anaconda3"):
                path = Path(os.path.expanduser(f"~/{candidate}/etc/profile.d/conda.sh"))
                if path.exists():
                    conda_sh = str(path)
                    break
            if not conda_sh:
                conda_sh = os.path.expanduser("~/miniconda3/etc/profile.d/conda.sh")
        gnina_activate = f"source {conda_sh} && conda activate {env_path}"

    ld_library_path = cfg.get("gnina_ld_library_path")
    if not ld_library_path and env_path:
        env_lib_paths: list[str] = []
        env_path_obj = Path(env_path)
        if env_path_obj.exists():
            torch_libs = sorted(
                env_path_obj.glob("lib/python*/site-packages/torch/lib")
            )
            env_lib_paths.extend(str(path) for path in torch_libs if path.is_dir())

            nvidia_libs = sorted(
                env_path_obj.glob("lib/python*/site-packages/nvidia/*/lib")
            )
            env_lib_paths.extend(str(path) for path in nvidia_libs if path.is_dir())

            lib_dir = env_path_obj / "lib"
            if lib_dir.is_dir():
                env_lib_paths.append(str(lib_dir))

        ld_library_path = _join_existing_library_paths(env_lib_paths)

    # Auto-detect from PyTorch or conda when no env_path is configured
    if not ld_library_path:
        ld_library_path = _auto_detect_cudnn_path()

    return gnina_activate, ld_library_path


def _join_existing_library_paths(paths: list[str]) -> str | None:
    """Join existing library directories into LD_LIBRARY_PATH string."""
    resolved: list[str] = []
    seen: set[str] = set()
    for raw in paths:
        path = Path(raw).expanduser()
        if not path.is_dir():
            continue
        normalized = str(path.resolve())
        if normalized in seen:
            continue
        seen.add(normalized)
        resolved.append(normalized)
    if not resolved:
        return None
    return os.pathsep.join(resolved)


def _auto_detect_cudnn_path() -> str | None:
    """Auto-detect LD_LIBRARY_PATH for GNINA from PyTorch or conda environments."""
    try:
        from hedgehog.setup._gnina import _collect_gnina_library_paths
    except ImportError:
        return None

    return _join_existing_library_paths(_collect_gnina_library_paths())


def _get_gnina_output_directory(cfg, base_folder):
    """Get GNINA output directory path."""
    gnina_config = cfg.get("gnina_config", {})
    cfg_out_dir = gnina_config.get("output_dir") or cfg.get("gnina_output_dir")
    if cfg_out_dir:
        out_dir_candidate = Path(cfg_out_dir)
        return (
            out_dir_candidate
            if out_dir_candidate.is_absolute()
            else (base_folder / out_dir_candidate)
        )
    return base_folder / "stages" / "05_docking" / "gnina"


def _resolve_receptor_path(receptor_pdb, base_folder=None):
    """Resolve receptor path to absolute, checking multiple locations.

    Args:
        receptor_pdb: Original receptor path from config
        base_folder: Base folder for relative path resolution

    Returns:
        Resolved Path object or None if not found
    """
    if not receptor_pdb:
        return None

    receptor_path = Path(receptor_pdb)
    if receptor_path.is_absolute() and receptor_path.exists():
        return receptor_path

    project_root = Path(__file__).parent.parent.parent.parent.parent

    if not receptor_path.is_absolute():
        candidate = (project_root / receptor_pdb).resolve()
        if candidate.exists():
            return candidate

    candidate = Path(receptor_pdb).resolve()
    if candidate.exists():
        return candidate

    return None
