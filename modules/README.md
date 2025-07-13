# Adding external modules

### Add MolScore module:
1. Clone `MolSCore` repo from source:
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
    export PYTHONPATH=$PYTHONPATH:/path/to/modules/MolGenBenchmark/modules/MolScore
    ```

### Add pyscreener module:
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

### Add MCE-18 module:
1. `MCE-18` implementation file was loaded from source: 
    [https://github.com/Tong-Du/MCE-18.git](https://github.com/Tong-Du/MCE-18.git)

    Download script and put it into `modules` folder.


When the evironment setup adjust `./configs` folder and specify configs based on metrics you want to calculate. And then run code:
```bash
python main.py
```


