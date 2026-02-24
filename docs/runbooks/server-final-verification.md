# Server Final Verification

Use this flow for fresh end-to-end validation on a target server.

```bash
ssh server
source ~/miniforge/etc/profile.d/conda.sh
conda activate base

cd ~/work/Projects
git clone https://github.com/LigandPro/hedgehog.git
cd hedgehog

uv sync

# --yes is not required (auto-accept by default)
uv run hedgehog setup aizynthfinder

# GNINA auto GPU is enabled by default
uv run hedgehog run --auto-install --out results/run_gpu_verify
```
