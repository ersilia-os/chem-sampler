# Generative Models Audit

Benchmarks Ersilia Model Hub generators on their own terms — output quality
(Tanimoto to seed, validity, timing) — independent of any annotator or
`hill_climb` run. Scripts run in numbered order; outputs land in `results/`,
model checkpoints in `models/`, per-model conda envs in `envs_cpu/` (all
gitignored). Script 03 needs `pip install -e ".[audit]"` (PyYAML) first.

## Scripts

| Script | Description |
|---|---|
| `scripts/01_prepare_chembl_splits.py` | Downloads ChEMBL, filters to single-component 250-450 Da compounds, samples 1000, and splits them into 100 ten-compound seed files (`results/ChEMBL_splits/split_*.csv`). |
| `scripts/02_fetch_generative_models.py` | Fetches the Ersilia Model Hub catalog, keeps the Ready generative (Sampling / Generation) models (Archived / In progress ones are skipped), and clones + `eosvc download`s each one's checkpoint into `--path-to-models` (default `models/`). |
| `scripts/03_build_cpu_envs.py` | Builds a CPU-only conda env per model (from its `install.yml` or `Dockerfile`) at `--path-to-envs/<model-id>-cpu` (default `envs_cpu/`). Builds with the user site off (`PYTHONNOUSERSITE=1`) so envs also work on machines with another home, e.g. the cluster. |
| `scripts/04_run_split.py` | Runs one model on one seed split inside its env (no conda needed) and writes `results/<model-id>/csv/split_NNN.csv` (`_repeated` with `--repeated`), plus a `.log` and a `.json` (wall time, host, CPUs, model commit) in `results/<model-id>/logs/`. Safe to re-submit: finished outputs are skipped. |
