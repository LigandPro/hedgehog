## Descriptors module

Molecular descriptors calculation and filtering for MolGenBenchmark. This stage computes per molecule physicochemical with RDKit.

The stage reads `configs/config_descriptors.yml` 

### What it does
- Computes per molecule metrics: chars of molecules, number of atoms, number of heavy atoms, number of hetero atoms, number of Nitrogen atoms, fraction of Nitrogen atoms, molecule charge (binary), molecular weight, logP, clogP (calculated), Sw, ring size, number of rings, number of aromatic rings, number of fused aromatic rings, number of rigid bonds, number of rotatable bonds, hydrogen bond donors (HBDs), hydrogen bond acceptors (HBAs), fraction sp3, topological polar surface area (TPSA), quantitative estimation of druglikeness (QED)
- Writes per molecule CSV, and pass and fail flags against configurable borders.


### Inputs
Supported upstream file formats (handled in preprocessing): CSV (any "smiles" column case), SMI/TXT (one SMILES per line, optional `,model_name`).