## Docking module

Prepares ligands and generates runnable scripts for Gnina and/or smina. 


### Config
Please extend config with proper information, fitting your task in 
`src/hedgehog/configs/config_docking.yml`:


**Run docking stage within entire pipeline:**\
Set `run: True` and adjust if needed [config_docking.yml](../../configs/config_docking.yml).

**Run docking stage only:**
```bash
uv run hedgehog run --stage docking

# Alternatively, using the short-name alias:
uv run hedge run --stage docking
```

### Running the generated scripts manually
Use prepared running scripts in `results` folder for your run:

```bash
./run_gnina.sh    
./run_smina.sh    
```

Notes:
- Scripts are executed with working directory set to `Docking`
- For Gnina, the output path is absolute by default to avoid nested path duplication.
- You can set `tools: gnina` or `tools: smina` to limit execution.
- Please take a look at `smina_example.ini` config to utilizr for `smina` docking. 
- 
- Set `auto_run: False` in the docking config to generate scripts without executing them.
- Use `gnina_extra_args` to pass flags like `--log score.log` or scoring options.


### Ligand Preparation
GNINA requires SDF format input files. The pipeline supports two modes for converting SMILES/CSV inputs to SDF:

#### With `ligand_preparation_tool`
If `ligand_preparation_tool` is configured in `src/hedgehog/configs/config.yml`, the pipeline uses the ligand preparation tool to convert input files:

- **Input format detection**: Automatically detects CSV (uses `-icsv`) or SMI files (uses `-ismi`)
- **Output**: Converts to SDF format saved to `docking/prepared_for_gnina.sdf`
- **Command**: `ligand_preparation_tool -icsv ligands.csv -osd docking/prepared_for_gnina.sdf`
- **Benefits**: Provides advanced preprocessing including stereo expansion, tautomerization, 3D geometry optimization, and other chemical transformations

To enable, add to `src/hedgehog/configs/config.yml`:
```yaml
ligand_preparation_tool: /path/to/ligand_preparation_tool
```

#### Without `ligand_preparation_tool`
If `ligand_preparation_tool` is **not** provided, the pipeline automatically falls back to built-in RDKit conversion:

1. Reads SMILES from the input CSV file
2. Converts SMILES to molecules using `Chem.MolFromSmiles()`
3. Adds hydrogens with `Chem.AddHs()`
4. Generates 3D conformations using `AllChem.EmbedMolecule()` with ETKDG method
5. Optimizes geometries with `AllChem.UFFOptimizeMolecule()`
6. Writes to `ligands_prepared.sdf`

### Common issues and fixes
- Gnina writes output into a duplicated folder path
  - Cause: giving Gnina a relative `-o` while launching from `Docking`. Gnina may prepend the working directory, resulting in duplicated segments.
  - Fix: we write an absolute output path (`gnina_output_dir` resolved to an absolute path). Ensure your config is up to date.
