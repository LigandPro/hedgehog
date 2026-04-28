# Adding external modules

## Add MCE-18 module:
1.  `MCE-18` implementation file was loaded from source: 
    [https://github.com/Tong-Du/MCE-18.git](https://github.com/Tong-Du/MCE-18.git)

    We thank authors for this brilliant open-source implementation of MCE-18 metric from https://pubs.acs.org/doi/10.1021/acs.jmedchem.9b00004 paper. 
<!-- 
## Install SYBA (required for synthesis stage)
SYBA is automatically installed when you create the conda environment from `environment.yml` (which includes the `lich` channel).

**If installing manually:**
```bash
conda activate hedge_env
mamba install lich::syba
```
Or with conda:
```bash
conda install -c lich syba
``` -->
<!-- 
**If conda installation fails**, you can install from source:
```bash
git clone https://github.com/lich-uct/syba.git
cd syba
conda activate hedge_env
pip install .
``` -->
<!-- Or uncomment the git pip dependency in `environment.yml` and recreate the environment.

## Install Eli Lilly Medchem Rules
Download Eli Lilly Medchem Rules via conda or mamba (preferred):
```bash
mamba install lilly-medchem-rules
```
Or with conda:
```bash
conda install -c conda-forge lilly-medchem-rules
``` -->



## Add AiZynthFinder retrosynthesis module 
### Automated Installation (Recommended)

**Quick setup via CLI:**
```bash
uv run hedgehog setup aizynthfinder
```

The setup command will automatically:
- Install the optional `retrosynthesis` dependency into the project environment
- Download public data 
- Set up logging configuration

Legacy fallback (script):

```bash
cd modules
./install_aizynthfinder.sh
```

### Manual Installation
If you prefer to install manually:

1. **Install the optional dependency**:
    ```bash
    uv sync --extra retrosynthesis
    ```

2. **Download public data**:
    ```bash
    mkdir -p modules/aizynthfinder/public modules/aizynthfinder/aizynthfinder/data
    uv run python -m aizynthfinder.tools.download_public_data modules/aizynthfinder/public
    cp src/hedgehog/synthesis/logging.yml modules/aizynthfinder/aizynthfinder/data/logging.yml
    ```
3. **Continue environment setup** following main [README.md](../README.md)

**Configure** in `src/hedgehog/configs/config_synthesis.yml`:
    - Set `run: True` to enable the stage
    - Adjust `nproc` for parallel processing if needed

The synthesis stage will run automatically after structural filters and before docking in the pipeline.
