import pandas as pd

from hedgehog.configs.logger import logger
from hedgehog.utils.datamol_import import import_datamol_quietly
from hedgehog.utils.input_paths import find_latest_input_source

dm = import_datamol_quietly()


def _find_latest_input_source(base_folder):
    """Find the most recent input source file for docking.

    Supports both new hierarchical structure and legacy flat structure.
    """
    path = find_latest_input_source(base_folder)
    if path:
        logger.debug("Using docking input: %s", path)
    return path


def _prepare_ligands_dataframe(df, output_csv):
    """Prepare ligands CSV from input DataFrame with SMILES validation."""
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    skipped_smiles = []

    for _, row in df.iterrows():
        smi = str(row["smiles"])
        try:
            mol = dm.to_mol(smi)
            if mol is None:
                skipped_smiles.append(smi)
                continue
        except Exception:
            skipped_smiles.append(smi)
            continue

        model_name = str(row["model_name"])
        mol_idx = str(row["mol_idx"])

        rows.append(
            {
                "smiles": smi,
                "name": mol_idx,
                "model_name": model_name,
                "mol_idx": mol_idx,
            }
        )

    output_df = pd.DataFrame(rows, columns=["smiles", "name", "model_name", "mol_idx"])
    output_df.to_csv(output_csv, index=False)

    if skipped_smiles:
        skip_path = output_csv.parent / "skipped_smiles.txt"
        with open(skip_path, "w") as f:
            for smi in skipped_smiles:
                f.write(f"{smi}\n")
        logger.warning(
            "Some SMILES could not be parsed for docking: %d/%d. See %s",
            len(skipped_smiles),
            len(df),
            skip_path,
        )
    return {
        "csv_path": str(output_csv),
        "total": len(df),
        "written": len(rows),
        "skipped": len(skipped_smiles),
    }
