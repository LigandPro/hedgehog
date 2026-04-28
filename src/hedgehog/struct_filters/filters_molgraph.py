"""Molecular graph statistics filter."""

import pandas as pd

from hedgehog._constants import CFG_STRUCT_FILTERS
from hedgehog.configs.logger import load_config, logger
from hedgehog.struct_filters._helpers import (
    _resolve_scheduler,
    add_model_name_col,
    mc,
)
from hedgehog.utils.datamol_import import import_datamol_quietly
from hedgehog.utils.parallel import resolve_n_jobs

dm = import_datamol_quietly()
_MOLGRAPH_CATALOG = None
_MOLGRAPH_SEVERITY_BY_ENTRY: list[int] | None = None


def _get_molgraph_catalog_and_severity():
    global _MOLGRAPH_CATALOG
    global _MOLGRAPH_SEVERITY_BY_ENTRY

    if _MOLGRAPH_CATALOG is None or _MOLGRAPH_SEVERITY_BY_ENTRY is None:
        graph_path = mc.utils.loader.get_data_path("graph.csv")
        graph_df = pd.read_csv(graph_path)
        _MOLGRAPH_SEVERITY_BY_ENTRY = graph_df["severity"].astype(int).tolist()
        _MOLGRAPH_CATALOG = mc.catalogs.NamedCatalogs.unstable_graph(
            severity_threshold=1
        )
    return _MOLGRAPH_CATALOG, _MOLGRAPH_SEVERITY_BY_ENTRY


def _molgraph_max_severity_for_mol(mol):
    catalog, severity_by_entry = _get_molgraph_catalog_and_severity()
    matched_entries = catalog.GetMatches(mol)
    max_severity = 0
    for entry in matched_entries:
        try:
            entry_idx = int(entry.GetDescription())
        except (TypeError, ValueError):
            continue
        if 0 <= entry_idx < len(severity_by_entry):
            max_severity = max(max_severity, int(severity_by_entry[entry_idx]))
    return max_severity


def _compute_molgraph_max_severities(mols, n_jobs, scheduler):
    if len(mols) == 0:
        return []
    return dm.parallelized(
        _molgraph_max_severity_for_mol,
        mols,
        n_jobs=n_jobs,
        scheduler=scheduler,
        progress=False,
    )


def apply_molgraph_stats(config, mols, smiles_model_name_mols=None, progress_cb=None):
    logger.info("Calculating Molecular Graph statistics...")
    total_mols = len(mols)
    config_sf = load_config(config[CFG_STRUCT_FILTERS])
    n_jobs = resolve_n_jobs(config_sf, config)
    scheduler = _resolve_scheduler(config_sf, "molgraph_scheduler")
    logger.info("MolGraph workers: %d", n_jobs)

    max_severities = _compute_molgraph_max_severities(mols, n_jobs, scheduler)
    results = {"mol": mols, "molgraph_max_severity": max_severities}
    if progress_cb is not None:
        progress_cb(total_mols, max(1, total_mols))

    for s in range(1, 12):
        results[f"pass_{s}"] = [sev < s for sev in max_severities]
    logger.info("MolGraph derived pass_1..pass_11 from one severity catalog pass")
    results = pd.DataFrame(results)

    if smiles_model_name_mols is not None:
        results = add_model_name_col(results, smiles_model_name_mols)
    return results
