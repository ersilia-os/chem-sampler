# Generative Models Audit

Benchmarks Ersilia Model Hub generators on their own terms — output quality
(Tanimoto to seed, validity, timing) — independent of any annotator or
`hill_climb` run. Scripts run in numbered order; outputs land in `results/`,
model checkpoints in `models/`, per-model conda envs in `envs_cpu/` (all
gitignored). Script 03 and the map in script 07 need `pip install -e ".[audit]"` (PyYAML, stylia) first.

## Scripts

| Script | Description |
|---|---|
| `scripts/01_prepare_chembl_splits.py` | Downloads ChEMBL, filters to single-component 250-450 Da compounds, samples 1000, and splits them into 100 ten-compound seed files (`results/ChEMBL_splits/split_*.csv`). |
| `scripts/02_fetch_generative_models.py` | Fetches the Ersilia Model Hub catalog, keeps the Ready generative (Sampling / Generation) models (Archived / In progress ones are skipped), and clones + `eosvc download`s each one's checkpoint into `--path-to-models` (default `models/`). |
| `scripts/03_build_cpu_envs.py` | Builds a CPU-only conda env per model (from its `install.yml` or `Dockerfile`) at `--path-to-envs/<model-id>-cpu` (default `envs_cpu/`). Builds with the user site off (`PYTHONNOUSERSITE=1`) so envs also work on machines with another home, e.g. the cluster. |
| `scripts/04_run_split.py` | Runs one model on one seed split inside its env (no conda needed) and writes `results/<model-id>/csv/split_NNN.csv` (`_repeated` with `--repeated`), plus a `.log` and a `.json` (wall time, host, CPUs, model commit) in `results/<model-id>/logs/`. Safe to re-submit: finished outputs are skipped. |
| `scripts/05_submit_slurm.py` | Submits a Slurm array that runs one model over all 100 splits (task N = split N, plus one extra task re-running a split as `_repeated`), on the lab nodes only (`spot_cpu` + lab nodelist; refuses others). Node names and the submit host come from `CHEMSAMPLER_LAB_CPU_NODES`, `CHEMSAMPLER_LAB_GPU_NODES` and `CHEMSAMPLER_SSH_HOST`, not from the code. Slurm `.out`/`.err` go to `results/<model-id>/out/`. `--dry-run` prints the batch script. |
| `scripts/06_evaluate.py` | Evaluates one model's outputs: % compounds with >= 1 null / >= 10 nulls / all slots null, max / P95 / mean Tanimoto to the input (mean ± std over compounds), retention of the input (same Murcko scaffold, same generic scaffold, contains the whole input, heavier; plus a breakdown at the most specific level kept), invalid, duplicate, input-echo and artifact counts, time per 100 compounds, and the comparison of a `_repeated` run. Writes `results/<model-id>/analysis.{json,md}`. Needs rdkit + numpy. |
| `scripts/07_plot_examples.py` | Draws `results/<model-id>/png/examples_grid.png` (and `png/chemical_space.png`, a map on the eos1klk t-SNE of the Ersilia reference library): one row per compound in the model's `examples/run_input.csv` (3 for most models) with the input first (INPUT 1, INPUT 2, ...) and its 5 most similar outputs (labelled `TS=<Tanimoto similarity to the input> (top-k)`) from the shipped `examples/run_output.csv`, or from `--outputs`. The map has one square scatter per example compound, side by side (seed black, outputs crimson, dot area by local density), titled `model | slug | Type`, with no legend, drawn with stylia (print format, article style). Needs rdkit, Pillow and numpy, plus stylia for the map (the `audit` extra). |
