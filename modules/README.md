# Adding external modules

## Add MolScore module:
1. Clone `MolSCore` from source:
    ```bash
    git clone https://github.com/MorganCThomas/MolScore.git
    ```
2. Install `MolScore` into environment:
    ```bash
    cd MolScore
    pip install .
    ```
3. Export path to `MolScore` module:
    ```bash
    export PYTHONPATH=$PYTHONPATH:/path/to/modules/hedge/modules/MolScore
    ```

## Add pyscreener module:
1. Clone `pyscreener` repo from source:
    ```bash
    git clone https://github.com/coleygroup/pyscreener.git
    ```
2. Install nessesary libraries: 
    ```bash
    conda install autodock-vina mgltools -c conda-forge -c bioconda -y
    ```
3. Add *prepare_receptor* naming:
    ```bash
    cd $CONDA_PREFIX/bin && ln -s prepare_receptor4.py prepare_receptor
    ```
4. Test `pyscreener` installation:
    ```bash
    pyscreener --config modules/pyscreener/integration-tests/configs/test_vina.ini --smoke-test
    pyscreener --smoke-test --screen-type vina --metadata-template '{"software": "vina"}'
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
    mv ../../src/hedge/stages/synthesis/logging.yml aizynthfinder/data/
    ```

3. **Configure** in `configs/config_synthesis.yml`:
    - Set `run: True` to enable the stage
    - Adjust `nproc` for parallel processing if needed

4. The synthesis stage will run automatically after structural filters and before docking in the pipeline.