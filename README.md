# 🦔 HEDGE
**Hierarchical Evaluation of Drug GEnerators**

![HEDGE Pipeline](data/imgs/pipeline_structure.png)


Comprehensive benchmark pipeline for evaluating generative models in molecular design.

### 5-Stage Filtering Pipeline:
1) Physicochemical Descriptors: 22 descriptors ([descriptors folder](src/hedge/stages/descriptors/))
2) Structural Filters: 6 structural criteria with ~2500 SMARTS patterns ([structural filters folder](structFilters/))
3) Synthesis evaluation ([synthesis folder](synthesis/))
4) Docking: able to calculate docking score with smina and/or GNINA docking tools ([docking folder](docking/))
5) Medicinal Chemists evaluation: calculate 4 medichinal chemists evaluation criteria ([medChem folder](medChem/))

## Setup & Run
**Hybrid approach:** Conda + UV

```bash
# Clone repository
git clone https://github.com/LigandPro/hedge.git
cd hedge

# Create conda environment for Lilly binaries (one-time setup)
conda env create -f environment.yml

# Every time you work
conda activate hedge_env

# Run full pipeline
uv run python main.py run

# Run specific stage
uv run python main.py run --stage descriptors

# Get help
uv run python main.py --help
```
 
### Adding external modules
To run this benchmark properly you need to install some extra packages. Go to [modules folder](modules/) and follow `README.md` inside.



### Configure your run
Edit config for each stage in [configs folder](configs/) based on metrics you want to calculate.

<!-- ## AIZythFinder retrosynthesis
To run AIZythFinder retrosynthesis clone code from source:
```bash
git clone https://github.com/MolecularAI/aizynthfinder.git
```


## docking score
--- -->
<!-- ## REINVENT4 fine-tune
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
--- -->
