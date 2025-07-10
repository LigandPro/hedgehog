# MolGenBenchmark
Benchmarking Generative Models for Molecular Design

## Metrics Calculation
Clone `MolGenBenchmark` repository repo for source:
   ```bash
   git clone https://github.com/LigandPro/MolGenBenchmark.git
   ```

### Environment preparation
1. Create environment with `environment.yml` file
2. Install SYBA using conda or mamba (preferred):
      ```bash
      mamba install lich::syba
      ```
      or clone SYBA from source:
      ```bash
      git clone https://github.com/lich-uct/syba.git && cd syba && conda activate YOUR_ENV && pip install .
      ```
3. Download Eli Lilly Medchem Rules via conda or mamba (preferred):
   ```bash
   mamba install lilly-medchem-rules
   ```
      
### Adding external modules
To run this benchmark propertly you need to install some extra packages. Go to [modules folder](projects/MolGenBenchmark/modules/README.md) and follow `README.md` inside.


When the evironment setup adjust `./configs` folder and specify configs based on metrics you want to calculate. And then run code:
```bash
python main.py
```


## AIZythFinder retrosynthesis
To run AIZythFinder retrosynthesis clone code from source:
```bash
git clone https://github.com/MolecularAI/aizynthfinder.git
```


## docking score
---
## REINVENT4 fine-tune
To fine-tune REINVENT4 follow these steps:
1. Clone REINVENT4 repository and setup environment:
    ```bash 
    git clone https://github.com/MolecularAI/REINVENT4.git
    ```
2. Configure transfer learning (aka fine-tuning)
   1. Adjust `configs/toml/transfer_learning.toml` following provided `configs/toml/README.md` instructions, 
   2. Set input model file for Mol2Mol generator as provided by authors `priors/reinvent.prior`.
   3. Set the following parameters:
     ```ini
     num_epochs = 1000
     save_every_n_epochs = 10
     batch_size = [adjust appropriately to reduce training time]
     ```
3) Train the model:
     Run the training using the modified configuration file. It takes approximetely 72 hours to train a model on ~750 samples with that setup. 
     
     Once trained, fine-tuned model can be used for downstream evaluation and benchmarking tasks.
---
