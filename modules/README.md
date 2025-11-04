# Adding external modules

## Add MolScore module:
1. Clone `MolSCore` from source:
    ```bash
    cd modules
    git clone https://github.com/MorganCThomas/MolScore.git
    ```

## Add pyscreener module:
1. Clone `pyscreener` repo from source:
    ```bash
    git clone https://github.com/coleygroup/pyscreener.git
    ```

## Add MCE-18 module:
1. `MCE-18` implementation file was loaded from source: 
    [https://github.com/Tong-Du/MCE-18.git](https://github.com/Tong-Du/MCE-18.git)

    Download script and put it into `modules` folder.
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
1. **Clone the retrosynthesis repository**:
    ```bash
    cd modules
    git clone git@github.com:LigandPro/retrosynthesis.git
    ```

2. **Download public data** (model files):
    ```bash
    cd modules/retrosynthesis/aizynthfinder
    mkdir -p public data
    # Using uv run
    uv run python -m aizynthfinder.tools.download_public_data ./public
    mv ../../../src/hedge/stages/synthesis/logging.yml data/
    ```

3. **Configure** in `configs/config_synthesis.yml`:
    - Set `run: True` to enable the stage
    - Adjust `nproc` for parallel processing if needed

4. The synthesis stage will run automatically after structural filters and before docking in the pipeline.