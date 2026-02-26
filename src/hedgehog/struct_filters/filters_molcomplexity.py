"""Molecular complexity filters."""

import contextlib
import io

import pandas as pd
from rdkit import Chem

from hedgehog._constants import CFG_STRUCT_FILTERS
from hedgehog.configs.logger import load_config, logger
from hedgehog.struct_filters._helpers import (
    _silence_worker_stdio,
    add_model_name_col,
    mc,
)
from hedgehog.utils.parallel import parallel_map, resolve_n_jobs

_COMPLEXITY_FILTERS_CACHE = {}
_MOLCOMPLEXITY_ALERT_NAMES: list[str] | None = None


def _get_complexity_filter(metric_name):
    """Get cached complexity filter instance for worker process."""
    complexity_filter = _COMPLEXITY_FILTERS_CACHE.get(metric_name)
    if complexity_filter is None:
        complexity_filter = mc.complexity.ComplexityFilter(
            complexity_metric=metric_name
        )
        _COMPLEXITY_FILTERS_CACHE[metric_name] = complexity_filter
    return complexity_filter


def _compute_molcomplexity_one(args):
    """Compute all complexity metrics for one molecule."""
    mol_idx, mol = args
    alert_names = _MOLCOMPLEXITY_ALERT_NAMES or []
    row = {"_mol_idx": mol_idx}

    if mol is None:
        row.update({f"pass_{name}": False for name in alert_names})
        row["pass"] = False
        row["pass_any"] = False
        return row

    try:
        prepared_mol = Chem.RemoveHs(mol)
    except Exception:
        prepared_mol = mol

    passed_all = True
    passed_any = False

    for name in alert_names:
        complexity_filter = _get_complexity_filter(name)
        if name == "smcm":
            with (
                contextlib.redirect_stdout(io.StringIO()),
                contextlib.redirect_stderr(io.StringIO()),
            ):
                try:
                    passed = bool(complexity_filter(prepared_mol))
                except Exception:
                    passed = False
        else:
            try:
                passed = bool(complexity_filter(prepared_mol))
            except Exception:
                passed = False

        row[f"pass_{name}"] = passed
        passed_all = passed_all and passed
        passed_any = passed_any or passed

    row["pass"] = passed_all
    row["pass_any"] = passed_any
    return row


def _init_molcomplexity_worker(alert_names: list[str]) -> None:
    global _MOLCOMPLEXITY_ALERT_NAMES
    _MOLCOMPLEXITY_ALERT_NAMES = alert_names


def _init_molcomplexity_worker_quiet(alert_names: list[str]) -> None:
    _silence_worker_stdio()
    _init_molcomplexity_worker(alert_names)


def apply_molcomplexity_filters(config, mols, smiles_model_name_mols=None):
    logger.info("Calculating Complexity filters...")
    config_path = config.get(CFG_STRUCT_FILTERS)
    config_sf = load_config(config_path) if config_path else {}
    n_jobs = resolve_n_jobs(config_sf, config)
    logger.info("MolComplexity workers: %d", n_jobs)

    alert_names = mc.complexity.ComplexityFilter.list_default_available_filters()
    items = [(i, mol) for i, mol in enumerate(mols)]
    rows = parallel_map(
        _compute_molcomplexity_one,
        items,
        n_jobs,
        initializer=(
            _init_molcomplexity_worker_quiet
            if n_jobs > 1 and len(items) > 1
            else _init_molcomplexity_worker
        ),
        initargs=(alert_names,),
    )
    for row in rows:
        mol_idx = row.pop("_mol_idx")
        row["mol"] = mols[mol_idx]
    final_result = pd.DataFrame(rows)

    if smiles_model_name_mols is not None:
        final_result = add_model_name_col(final_result, smiles_model_name_mols)
    return final_result
