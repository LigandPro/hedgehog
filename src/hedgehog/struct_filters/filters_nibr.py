"""NIBR (Novartis) structural filters."""

import pandas as pd

from hedgehog._constants import CFG_STRUCT_FILTERS
from hedgehog.configs.logger import load_config, logger
from hedgehog.struct_filters._helpers import (
    add_model_name_col,
    mc,
    _ensure_dataframe_length,
    _resolve_scheduler,
    _silence_worker_stdio,
    _split_indexed_mols,
)
from hedgehog.utils.parallel import parallel_map, resolve_n_jobs


def _process_nibr_chunk(args):
    """Run NIBR filter for one chunk with single-worker medchem call."""
    indexed_chunk, keep_details = args
    mol_indices = [item[0] for item in indexed_chunk]
    mol_chunk = [item[1] for item in indexed_chunk]

    nibr_filters = mc.structural.NIBRFilters()
    chunk_df = nibr_filters(
        mols=mol_chunk,
        n_jobs=1,
        scheduler="threads",
        keep_details=keep_details,
    )
    if len(chunk_df) != len(indexed_chunk):
        logger.warning(
            "NIBR chunk returned %d rows for %d molecules.",
            len(chunk_df),
            len(indexed_chunk),
        )
        template = chunk_df.iloc[-1].to_dict() if len(chunk_df) > 0 else None
        chunk_df = _ensure_dataframe_length(chunk_df, len(indexed_chunk), template)

    chunk_df = chunk_df.reset_index(drop=True)
    if "mol" in chunk_df.columns:
        chunk_df = chunk_df.drop(columns=["mol"])
    chunk_df["_mol_idx"] = mol_indices
    return chunk_df.to_dict("records")


def apply_nibr_filter(config, mols, smiles_model_name_mols=None):
    logger.info("Calculating NIBR filter...")
    config_sf = load_config(config[CFG_STRUCT_FILTERS])
    n_jobs = resolve_n_jobs(config_sf, config)
    scheduler = _resolve_scheduler(config_sf, "nibr_scheduler")
    logger.info("NIBR workers: %d", n_jobs)

    indexed_mols = list(enumerate(mols))
    n_workers = max(1, min(n_jobs, len(indexed_mols)))
    chunked_mols = _split_indexed_mols(indexed_mols, n_workers)
    chunk_payloads = [(chunk, True) for chunk in chunked_mols]
    chunk_results = parallel_map(
        _process_nibr_chunk,
        chunk_payloads,
        n_workers,
        chunksize=1,
        initializer=_silence_worker_stdio if n_workers > 1 else None,
    )

    rows = []
    for chunk_rows in chunk_results:
        rows.extend(chunk_rows)
    results = pd.DataFrame(rows)

    expected_indices = set(range(len(mols)))
    observed_indices = (
        set(results["_mol_idx"].tolist()) if "_mol_idx" in results else set()
    )
    if len(results) != len(mols) or observed_indices != expected_indices:
        logger.warning(
            "NIBR chunked processing mismatch (got %d/%d rows). Falling back to medchem native parallel call.",
            len(results),
            len(mols),
        )
        nibr_filters = mc.structural.NIBRFilters()
        results = nibr_filters(
            mols=mols,
            n_jobs=n_jobs,
            scheduler=scheduler,
            keep_details=True,
        )
    else:
        results = results.sort_values("_mol_idx").reset_index(drop=True)
        results["mol"] = [mols[i] for i in results["_mol_idx"]]
        results = results.drop(columns=["_mol_idx"])

    if smiles_model_name_mols is not None:
        results = add_model_name_col(results, smiles_model_name_mols)
    return results
