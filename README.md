# MolGenBenchmark
Benchmarking Generative Models for Molecular Design

## Metrics Calculation
Environment preparation
1. Go to `modules` folder and clone the following packages as external modules:
  ```bash
  cd modules
  ```

2. Install dependencies: \
   2.1 Install SYBA using mamba / conda (preferred):
      ```bash
      mamba install lich::syba
      ```
      or clone SYBA from source:
      ```bash
      git clone https://github.com/lich-uct/syba.git && cd syba && conda activate YOUR_ENV && pip install .
      ```
    2.2 Clone MolSCore repo from source:
      ```bash
      git clone https://github.com/MorganCThomas/MolScore.git
      ```
    2.3 Download MCE-18 implementation file from source: ```https://github.com/Tong-Du/MCE-18.git```
     
Metrics calculation:
1) Clone REPO repo for source:
   ```bash
   git clone https://github.com/REPO/REPO.git
   ```
2) Adjust config.py file 
2) Run:
   ```bash
   python main.py \
   --generated_mols_path path/to/.csv/.txt/.smi/files/with/molecules \
   --path_to_save path/to/save/data \
   --config config.yml
   ```
---
## AIZythFinder retrosynthesis
To run AIZythFinder retrosynthesis clone code from source:
  ```bash
  git clone https://github.com/MolecularAI/aizynthfinder.git
  ```

---
## _smina_ score
To evaluate _smina_ score, run
   ```bash
   conda install -c conda-forge smina
   ```
---
## REINVENT4 fine-tune
To fine-tune REINVENT4 follow these steps:
1) Clone REINVENT4 repository and setup environment:
  ```bash 
  git clone https://github.com/MolecularAI/REINVENT4.git
  ```
2) Configure transfer learning (aka fine-tuning)
   1. Adjust ```configs/toml/transfer_learning.toml``` following provided ```configs/toml/README.md``` instructions, 
   2. Set input model file for Mol2Mol generator as provided by authors ```priors/reinvent.prior```.
   3. Set the following parameters:
     ```ini
     num_epochs = 1000
     save_every_n_epochs = 10
     batch_size = [adjust appropriately to reduce training time]
     ```
3) Train the model:
     Run the training using the modified configuration file. One trained, this fine-tuned model can be used for downstream evaluation and benchmarking tasks.
---
