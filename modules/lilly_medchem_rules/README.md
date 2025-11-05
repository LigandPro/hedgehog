# Lilly MedChem Rules Binaries

This directory vendors the command-line utilities from the [Lilly-Medchem-Rules](https://github.com/IanAWatson/Lilly-Medchem-Rules) project.

The executables (`mc_first_pass`, `tsubstructure`, `iwdemerit`, `mc_summarise`) are required by the `medchem.structural.lilly_demerits` module. Shipping the compiled binaries here allows the pipeline to run under `uv` without a conda environment. The binaries were built on 2025-11-05 using the default GNU toolchain provided in this environment.

> Original project copyright (c) Eli Lilly and Company. See the upstream repository for licensing terms.
