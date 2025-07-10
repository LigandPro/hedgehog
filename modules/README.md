# MolScore config usage

To write a MolScore config make sure you have a proper script to run streamlit. Follow `modules/MolScore/molscore/gui/molscore_config` and check the script, it should be the following: 
```sh
config="$(dirname "$0")/config.py"
streamlit run $config "$@"
```

If yes, then you are ready to create a custom config:
1. Run a GUI on your pc using this command:
    ```bash
    modules/MolScore/molscore/gui/molscore_config
    ```
2. Create your config 
3. Run your config 
    ```python
    from molscore import MolScore

    ms = MolScore(model_name='WRITE_NAME',
                task_config='molscore/configs/JSONname.json',
                budget=10000
    )
    while not ms.finished:
        scores = ms.score(SMILES)
    ```

    


