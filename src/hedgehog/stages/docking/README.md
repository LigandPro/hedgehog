## Docking module

Prepares ligands and generates runnable scripts for Gnina and/or smina. 


### Config
Please extend config with proper information, fitting your task in 
`src/hedgehog/configs/config_docking.yml`:


**Run docking stage within entire pipeline:**\
Set `run: True` and adjust if needed [config_docking.yml](/src/hedgehog/configs/config_docking.yml).

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
GNINA requires SDF format input files. The pipeline supports multiple modes for handling ligand inputs:

#### Mode 1: Already SDF Format (No Conversion)
If you explicitly provide an SDF file via the `{tool_name}_ligands` config option, the pipeline uses it directly without any conversion:

**Configuration** (in `config_docking.yml`):
```yaml
smina_ligands: /path/to/ligands.sdf    # For SMINA
gnina_ligands: /path/to/ligands.sdf   # For GNINA
```

- **Supported formats**: `.sdf`, `.sdf.gz`, `.osd`, `.mol2`
- **No preprocessing**: File is used as-is
- **Fastest option**: No conversion overhead
- **Note**: This mode only applies when using `{tool_name}_ligands` config. The default pipeline flow (reading from previous stages) always converts to CSV first, then processes through Modes 2-4.

#### Mode 2: Direct RDKit Conversion (1:1 Mapping)
When `prepare_ligands: false` is set in `config_docking.yml`, the pipeline uses built-in RDKit conversion with 1:1 molecule mapping:

```yaml
prepare_ligands: false  # Skip advanced preprocessing
```
**Process:**
1. Reads SMILES from `ligands.csv` file (created from previous pipeline stages) or as input file if provided via --mols
2. Converts SMILES to molecules using `Chem.MolFromSmiles()`
3. Adds hydrogens with `Chem.AddHs()`
4. Generates 3D conformations using `AllChem.EmbedMolecule()` with ETKDG method
5. Optimizes geometries with `AllChem.UFFOptimizeMolecule()`
6. Writes to `_workdir/ligands_prepared.sdf`
7. **Maintains 1:1 mapping**: One input molecule → One output conformer

**Benefits**: Fast, maintains exact molecule count, no external dependencies

#### Mode 3: Advanced Preprocessing with Ligand Preparation Tool
When `prepare_ligands: true` (default) and `ligand_preparation_tool` is configured, the pipeline uses the external ligand preparation tool:

```yaml
prepare_ligands: true   # Default: true
```

**Configuration** (in `src/hedgehog/configs/config.yml`):
```yaml
ligand_preparation_tool: /path/to/ligand_preparation_tool
```
- **Note**: Tool should convert file to SDF format saved to `_workdir/prepared_for_{tool_name}.sdf`

#### Mode 4: Fallback to RDKit (when `prepare_ligands: true` but no tool)
If `prepare_ligands: true` but `ligand_preparation_tool` is **not** provided or not found, the pipeline automatically falls back to built-in RDKit conversion (same as Mode 2):


### Common issues and fixes
- Gnina writes output into a duplicated folder path
  - Cause: giving Gnina a relative `-o` while launching from `Docking`. Gnina may prepend the working directory, resulting in duplicated segments.
  - Fix: we write an absolute output path (`gnina_output_dir` resolved to an absolute path). Ensure your config is up to date.