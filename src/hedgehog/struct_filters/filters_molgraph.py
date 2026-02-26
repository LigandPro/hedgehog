"""Molecular graph statistics filter."""

import pandas as pd

from hedgehog._constants import CFG_STRUCT_FILTERS
from hedgehog.configs.logger import load_config, logger
from hedgehog.struct_filters._helpers import (
    _resolve_scheduler,
    add_model_name_col,
    mc,
)
from hedgehog.utils.parallel import resolve_n_jobs


def apply_molgraph_stats(config, mols, smiles_model_name_mols=None):
    logger.info("Calculating Molecular Graph statistics...")
    severities = list(range(1, 12))
    config_sf = load_config(config[CFG_STRUCT_FILTERS])
    n_jobs = resolve_n_jobs(config_sf, config)
    scheduler = _resolve_scheduler(config_sf, "molgraph_scheduler")
    logger.info("MolGraph workers: %d", n_jobs)

    results = {"mol": mols}
    for s in severities:
        out = mc.functional.molecular_graph_filter(
            mols=mols,
            max_severity=s,
            n_jobs=n_jobs,
            scheduler=scheduler,
            progress=False,
            return_idx=False,
        )
        results[f"pass_{s}"] = out
    results = pd.DataFrame(results)

    if smiles_model_name_mols is not None:
        results = add_model_name_col(results, smiles_model_name_mols)
    return results
