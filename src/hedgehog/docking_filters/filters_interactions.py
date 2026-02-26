"""Protein-ligand interaction filter using ProLIF."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import Chem

from hedgehog.configs.logger import logger
from hedgehog.utils.parallel import resolve_n_jobs


def apply_interaction_filter(
    mols: list[Chem.Mol],
    protein_pdb: str | Path,
    config: dict[str, Any],
) -> pd.DataFrame:
    """
    Apply interaction filter using ProLIF.

    Checks for required protein-ligand interactions.

    Args:
        mols: List of RDKit molecules with 3D coordinates
        protein_pdb: Path to protein PDB file
        config: Filter configuration dict

    Returns:
        DataFrame with columns: mol_idx, n_hbonds, interactions, pass
    """
    import MDAnalysis as mda
    import prolif as plf

    min_hbonds = config.get("min_hbonds", 0)
    required_residues = config.get("required_residues", [])
    forbidden_residues = config.get("forbidden_residues", [])
    interaction_types = config.get(
        "interaction_types", ["HBDonor", "HBAcceptor", "Hydrophobic", "VdWContact"]
    )

    logger.info("Running interaction filter (min_hbonds=%d)", min_hbonds)

    # Load protein
    u = mda.Universe(str(protein_pdb))
    protein = plf.Molecule.from_mda(u)

    # Convert RDKit mols to ProLIF molecules
    plf_mols = []
    for mol in mols:
        mol_h = Chem.AddHs(mol, addCoords=True)
        plf_mols.append(plf.Molecule.from_rdkit(mol_h))

    # Calculate interactions (exclude FaceToFace due to edge case errors)
    safe_interactions = [it for it in interaction_types if it != "FaceToFace"]
    fp = plf.Fingerprint(interactions=safe_interactions)
    n_jobs = resolve_n_jobs(config)
    fp.run_from_iterable(plf_mols, protein, n_jobs=n_jobs)

    df_interactions = fp.to_dataframe(drop_empty=False)

    # Build results DataFrame
    results = []
    for i in range(len(mols)):
        if i < len(df_interactions):
            row = df_interactions.iloc[i]
            # Count H-bonds
            hbond_cols = [
                c for c in row.index if "HBDonor" in str(c) or "HBAcceptor" in str(c)
            ]
            n_hbonds = sum(row[hbond_cols].values) if hbond_cols else 0

            # Check required residues
            has_required = True
            for res in required_residues:
                res_cols = [c for c in row.index if res in str(c)]
                if res_cols and not any(row[res_cols].values):
                    has_required = False
                    break

            # Check forbidden residues
            has_forbidden = False
            for res in forbidden_residues:
                res_cols = [c for c in row.index if res in str(c)]
                if res_cols and any(row[res_cols].values):
                    has_forbidden = True
                    break

            passed = n_hbonds >= min_hbonds and has_required and not has_forbidden
            interactions_str = ",".join([str(c) for c in row.index if row[c]])
        else:
            n_hbonds = 0
            passed = min_hbonds == 0 and not required_residues
            interactions_str = ""

        results.append(
            {
                "mol_idx": i,
                "n_hbonds": n_hbonds,
                "interactions": interactions_str,
                "pass_interactions": passed,
            }
        )

    df = pd.DataFrame(results)
    logger.info(
        "Interaction filter: %d/%d passed",
        int(df["pass_interactions"].sum()),
        len(df),
    )
    return df
