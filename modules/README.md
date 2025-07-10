### Adding external modules
1. Clone MolSCore repo from source:
    ```bash
    git clone https://github.com/MorganCThomas/MolScore.git
    ```
    and export path to MolScore module:
    ```bash
    export PYTHONPATH=$PYTHONPATH:/path/to/modules/MolGenBenchmark/modules/MolScore
    ```



2. MCE-18 implementation file was loaded from source: 
    [https://github.com/Tong-Du/MCE-18.git](https://github.com/Tong-Du/MCE-18.git)


When the evironment setup adjust `./configs` folder and specify configs based on metrics you want to calculate. And then run code:
```bash
python main.py
```


