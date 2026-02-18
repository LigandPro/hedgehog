## Descriptors Module

Computes 22 physicochemical descriptors per molecule using RDKit. 

The stage reads [config_descriptors](../../configs/config_descriptors.yml) file.

Descriptors:
- **Molecular properties**: chars, number of atoms, heavy atoms, hetero atoms.
- **Chemical features**: molecule charge (binary), molecular weight, logP, clogP, Sw.
- **Structure**: number of Nitrogen atoms, fraction of Nitrogen atoms, ring size, number of rings, number of aromatic rings, number of fused aromatic rings, number of rigid bonds, number of rotatable bonds, hydrogen bond donors (HBDs), hydrogen bond acceptors (HBAs), fraction sp3, topological polar surface area (TPSA).
- **Drug-likeness**: quantitative estimation of druglikeness (QED).

### Usage
Physicochemical Descriptors are calculated **twice** in **🦔 HEDGEHOG**: as a first stage (*could be tuned off*) and after any last stage of a pipeline to collect final molecules statistics (*counld not be turned off*).

**Run descriptors stage within entire pipeline:**\
Set `run: True` and adjust if needed [config_descriptors.yml](../../configs/config_descriptors.yml).

**Run descriptors stage only:**
```bash
uv run hedgehog run --stage descriptors

# Alternatively, using the short-name alias:
uv run hedge run --stage descriptors
```

**Note:** Ensure [config.yml](../../configs/config.yml) and [config_descriptors.yml](../../configs/config_descriptors.yml) are properly configured.
