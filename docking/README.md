## Docking module

Prepares ligands from latest pass list and generates runnable scripts for Gnina and/or smina using Pyscreener. 

### Inputs (auto-detected priority)
1) `Retrosynthesis/passRetrosynthesisSMILES.csv`
2) `StructFilters/passStructFiltersSMILES.csv`
3) `Descriptors/passDescriptorsSMILES.csv`
4) `sampledMols.csv`

### Config
Please extend config with proper information, fitting your task in 
`configs/config_docking.yml`:


### How to use in pipeline
1) Ensure updated config
2) Set `run: True` in `configs/config_docking.yml` and correct all absolute paths.
3) Run the main pipeline:

### Running the generated scripts manually
From the `Docking` directory:

```bash
./run_gnina.sh    # writes poses to gnina_results/gnina_out.sdf
./run_smina.sh    # launches pyscreener
```

Notes:
- Scripts are executed with working directory set to `Docking`
- For Gnina, the output path is absolute by default to avoid nested path duplication.
- You can set `tools: gnina` or `tools: smina` to limit execution.
- Set `auto_run: False` in the docking config to generate scripts without executing them.
- Use `gnina_extra_args` to pass flags like `--log score.log` or scoring options.


### Common issues and fixes
- Gnina writes output into a duplicated folder path
  - Cause: giving Gnina a relative `-o` while launching from `Docking`. Gnina may prepend the working directory, resulting in duplicated segments.
  - Fix: we write an absolute output path (`gnina_output_dir` resolved to an absolute path). Ensure your config is up to date.
