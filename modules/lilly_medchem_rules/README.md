# Lilly MedChem Rules binaries

Vendored CLI tools from
[Lilly-Medchem-Rules](https://github.com/IanAWatson/Lilly-Medchem-Rules)
so Hedgehog can run Lilly demerit scoring under `uv` without a conda env.

**Binaries here:** `mc_first_pass`, `tsubstructure`, `iwdemerit`, `mc_summarise`  
Used by `medchem.structural.lilly_demerits` when Stage 3 has:

```yaml
calculate_lilly: true
filter_lilly: true          # hard gate (optional)
lilly_demerit_cutoff: 160
```

Hedgehog prepends `modules/lilly_medchem_rules/bin` to `PATH` when needed.

Built 2025-11-05 with the default GNU toolchain in this environment.

> Original copyright (c) Eli Lilly and Company — see the upstream repo for license terms.
