import json
import uuid
from datetime import datetime

from hedgehog._constants import TOOL_GNINA, TOOL_MATCHA, TOOL_SMINA
from hedgehog.configs.logger import logger


def _generate_job_id(tool="dock"):
    """Generate a unique job ID."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    unique_id = uuid.uuid4().hex[:8]
    return f"{tool}_{timestamp}_{unique_id}"


def _save_job_metadata(
    ligands_dir,
    source_file,
    num_ligands,
    receptor_pdb,
    tools_prepared,
    scripts_prepared,
    ligands_csv,
    ligands_stats,
    job_ids,
    overall_job_id,
):
    """Save job metadata to JSON file."""
    ligands_dir.mkdir(parents=True, exist_ok=True)
    metadata = {
        "job_id": overall_job_id,
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "source_file": str(source_file),
        "num_ligands": num_ligands,
        "receptor_pdb": receptor_pdb,
        "tools_prepared": tools_prepared,
        "scripts": scripts_prepared,
        "ligands_csv": str(ligands_csv),
        "ligands_counts": ligands_stats,
        "jobs": {
            tool: {
                "name": tool,
                "job_id": job_id,
                "script": str(ligands_dir / "_workdir" / f"run_{tool}.sh"),
            }
            for tool, job_id in job_ids.items()
        },
    }
    meta_path = ligands_dir / "job_meta.json"
    with open(meta_path, "w") as f:
        json.dump(metadata, f, indent=2)


def _save_job_ids(ligands_dir, overall_job_id, job_ids):
    """Save job IDs to a simple text file."""
    ids_path = ligands_dir / "job_ids.txt"
    try:
        with open(ids_path, "w") as f:
            f.write(f"overall: {overall_job_id}\n")
            for tool in (TOOL_SMINA, TOOL_GNINA, TOOL_MATCHA):
                f.write(f"{tool}: {job_ids.get(tool, '')}\n")
    except Exception as e:
        logger.warning("Failed to write job_ids.txt: %s", e)


def _update_metadata_with_run_status(ligands_dir, run_status):
    """Update job metadata with run status."""
    meta_path = ligands_dir / "job_meta.json"
    try:
        metadata = {}
        if meta_path.exists():
            with open(meta_path) as f:
                metadata = json.load(f)
        metadata["run_status"] = run_status
        with open(meta_path, "w") as f:
            json.dump(metadata, f, indent=2)
    except Exception as e:
        logger.warning("Failed to update metadata with run status: %s", e)


def _parse_tools_config(cfg):
    """Parse and validate an explicit docking engine selection."""
    tools_cfg = cfg.get("tools", [TOOL_SMINA, TOOL_GNINA])

    if isinstance(tools_cfg, str):
        requested = [item.strip().lower() for item in tools_cfg.split(",")]
    elif isinstance(tools_cfg, (list, tuple)):
        requested = [str(item).strip().lower() for item in tools_cfg]
    else:
        raise ValueError("docking tools must be a string or a list")

    requested = [item for item in requested if item]
    if not requested:
        raise ValueError("docking tools must select at least one engine")

    aliases = {
        "all": (TOOL_SMINA, TOOL_GNINA, TOOL_MATCHA),
        "both": (TOOL_SMINA, TOOL_GNINA),
    }
    supported = {TOOL_SMINA, TOOL_GNINA, TOOL_MATCHA}
    selected: list[str] = []
    for item in requested:
        expanded = aliases.get(item, (item,))
        for tool in expanded:
            if tool not in supported:
                raise ValueError(f"Unsupported docking tool: {tool}")
            if tool not in selected:
                selected.append(tool)
    return selected
