---
name: gen-model-audit
description: >
  Audit one Ersilia Model Hub generative model (Task Sampling, Subtask Generation) end to end on the
  shared 1000-compound ChEMBL benchmark: set up its environment, read its code and summarise what it
  does, run it as a Slurm array on the lab nodes (100 splits of 10 compounds plus one repeated split),
  evaluate nulls, Tanimoto, retention of the input and output quality, draw the example grid and the
  chemical-space map, optionally measure a GPU, and fold the result into the status table of
  ersilia-os/ersilia issue #1919. Never pushes, tags or releases a model, and posts to GitHub only
  after the user has seen the text. Triggers include: "audit a generative model", "gen-model-audit",
  "/gen-model-audit", "benchmark eosXXXX", "run the 1000-compound benchmark on", "update the gen
  models table", "refresh issue 1919".
---

# Generative model audit

One model per run. Everything lives in `examples/gen-models-audit/` of `ersilia-os/chem-sampler` (public repo, branch `chemsampler-update` while it is being built). Results and scratch are gitignored.

## The rules (never break)

1. **Never push, tag or release** anything in an `ersilia-os/<model>` repo, and never change a version, without the user's explicit say-so for that model. Report problems and proposals instead. Releases are the maintainers' decision (output-visible changes are MAJOR in the model template).
2. **Lab nodes only.** Other nodes are billed to the PI. Script 05 refuses nodes outside `CHEMSAMPLER_LAB_CPU_NODES` / `CHEMSAMPLER_LAB_GPU_NODES`; never pass `--allow-any-node` on your own. GPU jobs use the partition and lab GPU nodes with `--nodes=1 --ntasks=1 --gres=gpu:1`, the same header pattern as the user's other array scripts.
3. **Do not poll Slurm.** After submitting, give the job id and the two monitoring commands, then stop; the user says when the jobs are done.
4. **`PYTHONNOUSERSITE=1` for every environment build and run.** Environments must not lean on `~/.local`. There is no conda on the cluster nodes: call the env's interpreter directly.
5. **No cluster address or node name in committed code or docs** (the repo is public). They live in environment variables: `CHEMSAMPLER_SSH_HOST`, `CHEMSAMPLER_LAB_CPU_NODES`, `CHEMSAMPLER_LAB_GPU_NODES` (set in the user's `~/.bashrc`; export them inline in a Claude session if the shell predates them).
6. **Ask before acting.** The repo's CLAUDE.md: plan mode for anything non-trivial, `AskUserQuestion` for choices, session notes in `.claude/` at the end. Scratch work goes in the repo's gitignored `tmp/<topic>/`.
7. **Post to GitHub only after the user has seen the final text.** `gh` is a snap and cannot read files on the shared lab filesystem: pipe the body (`cat final.md | gh issue edit 1919 -R ersilia-os/ersilia --body-file -`) and verify the live body afterwards. Never touch other people's comments.
8. **Never write a SMILES from memory**, and never type a number into a report that a file can give you.
9. **Keep table columns** unless the user names the column to drop.
10. Commit or push in chem-sampler only when asked, staging by file name.

## What the audit established (state at the end of the pilot)

- **Benchmark:** 1000 ChEMBL 37 compounds (single component, 250-450 Da, random sample, seed 42), 100 files of 10 (`results/ChEMBL_splits/split_001..100.csv`, column `smiles`). One Slurm task per split plus one extra task re-running `split_001` as `split_001_repeated.csv`.
- **Definitions:** duplicate = same canonical *isomeric* SMILES within a compound; input echo = output equal to its input ignoring stereochemistry (must never happen); null = empty output slot; Tanimoto = Morgan radius 2, 2048 bits, between each valid output and its own input, per compound max / P95 / mean, then mean ± sample std over compounds; retention = share of a compound's outputs keeping the input's Murcko scaffold / generic scaffold / the whole input, plus a breakdown at the most specific level kept.
- **Upstream fixes already pushed (unreleased):** eos9taz, eos633t, eos8fma, eos6ost, eos69e6 (1000 outputs before) and eos57bx (500) now return at most 100 molecules with no duplicates or input echo; eos4q1a and eos9p57 no longer return duplicates (K-Means picked the same molecule for two centres). 14 of the 19 models have output-affecting commits after their latest release (`Critical change after release` in `gen-models-master.csv`).
- **Type scheme (7 classes):** Fragment-seeded growth, Fragment recombination, Scaffold-based, Direct edit, Similarity-conditioned translation, Descriptor-conditioned, Reaction-based. Per model in `gen-models-master.csv`.
- **Metadata rule:** every generative model must have Output Consistency = Variable. Still Fixed: eos9taz, eos57bx, eos694w, eos5j3l, eos935d.
- **Pilot eos2401 (reference run):** 30.1% of compounds have at least one null, 14.6% at least ten, 10.3% return nothing; max / P95 / mean Tanimoto 0.22 / 0.18 / 0.13; 99.9% of outputs keep none of the whole input, Murcko or generic scaffold; fails on all single-ring-system and acyclic inputs; 21.5 min per 100 compounds on 4 CPUs; a second run shares only 1.6 molecules per compound but has the same statistics.
- **GPU (eos2401 only):** transformer step 15x (RTX 3090) and 23x (H200) faster than the CPU node; 1000 compounds in about 13 min on one 3090 against 3.6 h on the CPU array; output statistically equivalent. The library's joblib pool start-up costs a one-off 11-16 s per process. Details in `results/eos2401/md/gpu_assessment.md`.
- **Issue #1919** was rewritten and posted on the basis of this: wide status table, every model other than eos2401 "pending". The table is regenerated, not edited by hand (see step 11).

## Repository map

| Script | Does | Key flags |
|---|---|---|
| `01_prepare_chembl_splits.py` | ChEMBL download, filter, seed-42 sample, 100 split files. **No argparse: running it with `--help` runs the whole download. Run once, never to "test".** | none |
| `02_fetch_generative_models.py` | Catalog filter (Ready generative), clone/pull each repo, `eosvc download` | `--models id...` (overrides the catalog), `--path-to-models` |
| `03_build_cpu_envs.py` | Per-model conda env `envs_cpu/<id>-cpu` from `install.yml` or Dockerfile, user site off, warns on missing requirements | `--models`, `--force` |
| `04_run_split.py` | One model on one split with the env's interpreter; csv + log + provenance json; skips finished outputs; exit 1 if rows != inputs | `--model`, `--split-id`, `--repeated`, `--cpus`, `--force` |
| `05_submit_slurm.py` | Slurm array on the lab nodes (task N = split N, last task = repeated split) | `--model`, `--array`, `--time`, `--mem`, `--cpus`, `--max-parallel`, `--dry-run` |
| `06_evaluate.py` | Nulls, Tanimoto, retention, quality counts, time per 100 compounds, repeated-split comparison | `--model`, `--null-threshold` |
| `07_plot_examples.py` | `png/examples_grid.png` (3 x 6 molecules: input + top 5 outputs labelled `TS=<Tanimoto> (top-k)`, no legend, scaffold red / fragment blue by Type; the fragment rule is chosen by hand per model: by default the fragments kept by any drawn output, and for models whose outputs all grow from one fragment of the input, listed in `SHARED_FRAGMENT_MODELS` in the script (eos8zvb), only that parent ring system) and `png/chemical_space.png` (eos1klk t-SNE, one square scatter per seed side by side: seed black, outputs crimson, dot area by density, centred `model | slug | Type` title, no legend, stylia print / article) | `--model`, `--plots {both,grid,space}`, `--outputs`, `--highlight` |

Result layout per model: `results/<id>/csv/` (split outputs), `logs/` (model logs, provenance `.json`, submitted `.sbatch`), `out/` (Slurm `.out`/`.err`), `md/` (`analysis.md`, `summary.md`, `gpu_assessment.md`), `png/`, and `analysis.json` at the top.
Scratch generators in `tmp/`: `issue-1919/build_draft.py` (issue body from files), `eos2401-gpu/` (CPU phase profile harness, GPU env build, job files: eos2401-specific templates).

## One-time setup

1. `pip install -e ".[audit]"` in the chemsampler conda env: PyYAML for script 03 and `stylia` (pinned in the extra) for the map in script 07; rdkit, numpy, matplotlib and Pillow come with the core dependencies. Without `stylia`, `07 --plots grid` still works but `--plots space` fails on import.
2. Environment variables above in `~/.bashrc`.
3. Seeds: `results/ChEMBL_splits/` exists (script 01 is done; do not rerun).
4. For the chemical-space plot: `02 --models eos1klk`, `03 --models eos1klk`; the first `07 --plots space` projects 20,000 reference-library compounds (lazy-chemvis `smiles_100k.csv`, seed 42) and caches them in `cache/eos1klk/` (under a minute). Only eos1klk's t-SNE surrogate is run (Morgan fingerprint radius 2, 2048 bits, into one XGBoost model): the coordinates are identical to the full model's, which also computes RDKit descriptors, PCA, TMAP and UMAP and is about 300x slower. Do not run the full eos1klk `run.sh` for this.

## Per-model procedure

Do the steps in order; each ends in a stop or a check. Replace `<id>`.

1. **Fetch and build.** `python scripts/02_fetch_generative_models.py --models <id>` then `python scripts/03_build_cpu_envs.py --models <id>`. Check: no WARNING lines left, `PYTHONNOUSERSITE=1 envs_cpu/<id>-cpu/bin/python -m pip check` clean, and the model's own `examples/run_input.csv` runs through `run.sh` with that env first on PATH. `eosvc` exiting non-zero on an empty `model/framework/fit` is normal, not a failed checkpoint. Models with a non-conda step (eos4qda: opam, needs sudo) need the user.
2. **Read the code and write `results/<id>/md/summary.md`**: what the model does, in short bullets (input handling, how it generates, what it filters, output budget, failure modes, number of outputs). Check what it downloads at run time (`from_pretrained`, `hf_hub`, `urlretrieve`, `requests`): if it needs network, copy the cache to `cache/huggingface` (script 04 then sets `HF_HOME` and `HF_HUB_OFFLINE=1`).
3. **Local smoke test and sizing.** `python scripts/04_run_split.py --model <id> --split-id 1` (for a slow model, a 2-compound split in a tmp dir via `--path-to-splits`). Note wall time and memory, and which compounds give empty rows. Choose `--time`, `--mem`, `--cpus` for step 5. Estimates per 10-compound task: eos2401 about 2 min; eos6a1h about 50 min; eos69e6 about 35 min; eos55vx about 3-4 h (about 20 min per compound). Everything else untimed.
4. **Cluster smoke test.** `python scripts/05_submit_slurm.py --model <id> --array 1 --time ... --mem ...` (use `--dry-run` first). One task; check `results/<id>/csv/split_001.csv`, the `.json` (exit 0, host is a lab node) and an empty `.err`. A first task that fails on imports means the env leaked the user site: fix the env, not the job.
5. **Full array.** `python scripts/05_submit_slurm.py --model <id> --time ... --mem ...` (array `1-101%50`; a finished split is skipped, so the same command resumes after a failure; to change the limit of an array already submitted: `ssh $CHEMSAMPLER_SSH_HOST scontrol update JobId=<id> ArrayTaskThrottle=<n>`). **Give the job id, `ssh $CHEMSAMPLER_SSH_HOST squeue -j <id>` and `sacct -j <id> -X`, then stop.**
6. **When the user says the jobs are done, verify:** 100 `split_NNN.csv` plus `split_001_repeated.csv`, no `*.part.csv`, no non-zero `exit_code` in `logs/*.json`, empty `out/*.err`. Re-run step 5 for missing splits. Count rows with the csv module, not `wc -l` (files may lack a trailing newline).
7. **Evaluate.** `python scripts/06_evaluate.py --model <id>`; re-run whenever more splits arrive (it rewrites `analysis.json` and `md/analysis.md`).
8. **Plots.** `python scripts/07_plot_examples.py --model <id>`; view both images and check blanks, highlights and legend. Use `--outputs` to draw a fresh run instead of the shipped `examples/run_output.csv`.
9. **Quality review**, using 06's counts and the repo: invalid, duplicate, echo, multi-component, atom-map, dummy-atom must be 0 (stereoisomer-only repeats are allowed under the duplicate definition; report them); `run_output.csv` shape vs `run_columns.csv` and metadata Output Dimension; Description / Interpretation / Tag vs the measured behaviour; Output Consistency must be Variable; output size vs description; `gh run list -R ersilia-os/<id>` for red workflows; commits after the latest release (`gh api repos/ersilia-os/<id>/compare/<tag>...main`). List problems with evidence; propose fixes; do not push them.
10. **Optional: CPU phase profile and GPU test.** Only when the model has a neural sampling step. Time each phase on one split on a CPU node first (wrap the model's calls with timers, two passes in one process, watch for one-off pool or model-load costs). If the GPU-capable step is more than about half of the time, build a GPU env (the CPU env's `pip freeze` without torch, then the CUDA wheel of the same torch version; the lab GPU nodes' driver supports CUDA 12.8) and compare CPU, 3090 and H200 on `split_001` plus a 100-compound equivalence run checked with script 06. `tmp/eos2401-gpu/` has the harness, job templates and `build_assessment.py`. Report the model-side change as a proposal.
11. **Fold into the status table.** Update the model's row in `gen-models-master.csv` if its release, flag or size changed, run `tmp/issue-1919/build_draft.py` (it reads `analysis.json` for every benchmarked model and writes `tmp/issue-1919/draft.md`; unbenchmarked models stay "pending"), show the user the new text, and post it only after approval (rule 7).
12. **Session notes.** Append what changed, the real numbers, mistakes and open items to the dated `.claude/*-session-notes.md`.

### Refreshing the clones (any time, between models)

`python scripts/02_fetch_generative_models.py` with no `--models` pulls all Ready generative repos and re-runs `eosvc download`. Record `git rev-parse HEAD` per clone before, compare after, and read the new commits per model: `ersilia-bot` "updating metadata/readme [skip ci]" commits are CI follow-ups, anything else is a real change (code, columns, examples). Compare the catalog (`https://catalog.ersilia.io/api/models`: Release, Output Dimension, Title, Slug, Status) with `gen-models-master.csv` and flip "Critical change after release" for any model whose code changed after its latest release. The script's "N/19 checkpoint downloads confirmed" counts `eosvc` errors on empty folders too: check the checkpoint folders have content instead of trusting it.

## Decisions that belong to the user (ask, don't assume)

Releases and version bumps; any push to a model repo (code fixes, metadata, examples); changing a model's Output Consistency or output size; posting or editing the issue; committing or pushing in chem-sampler; GPU image or device changes in a model; running a very slow model in full; adding a dependency.

## Pitfalls already paid for

- **User-site leakage:** pip skipped packages found in `~/.local`, so 10 of 19 envs lacked jinja2, flask, sqlalchemy, ...; they failed on cluster nodes. Script 03 now sets `PYTHONNOUSERSITE=1` and warns on `pip check`.
- **`conda run` buffers output and does not exist on the nodes:** use `<env>/bin/python` (scripts 04 and 05 do).
- **`--nodelist` with `--nodes=1`** (several lab nodes listed) works as "any of these" on this cluster; the user's own array scripts use it.
- **`gh` snap cannot read paths on the shared lab filesystem:** `--body-file -` with stdin.
- **`pkill -f` with a pattern that appears in your own command line kills your shell:** kill by PID.
- **`--help` on a script without argparse runs the script** (script 01).
- **RDKit `MakeScaffoldGeneric` raises on hypervalent atoms** (P, S become 5-6-valent carbons): 06 treats that scaffold as having no generic form.
- **Examples files may lack a trailing newline:** count rows with csv.
- **A model's example inputs can include compounds it cannot handle** (empty rows are data, not errors).
- **Safe/joblib pool start-up** (eos2401): the first call that decodes more than 100 sequences, or the substructure filter, starts a process pool (11 s on the CPU node, 16 s on a GPU node); serial replacements remove it.
- **`analysis.json` goes stale** when splits are added; re-run 06 before quoting numbers.
- **Disk:** the shared home filesystem is about 99% full; a CUDA env is 11 GB. Reuse caches, do not copy envs.
- **A model can ignore its input without raising anything** (eos8zvb: an `AttributeError` swallowed by `except Exception: continue` sent every input to a benzene fallback): in the quality review compare each output's Tanimoto to its own input with its Tanimoto to a shuffled other input (`tmp/eos8zvb/input_dependence.py`); equal values mean the output does not depend on the input. Also call the model's own helper on a few trivial inputs and look at what it returns, and list which fragments / steps raise inside a bare `except` (eos8zvb had a second one: a `KekulizeException` made it skip the larger fragment, e.g. a purine, of a third of the inputs).
- **SMILES-unique is not diverse** (eos8zvb: the growth generator yields every intermediate stage, so about 80% of the 100 outputs are a substructure of another output, about 20 independent molecules per compound): in the quality review also report the share of outputs contained in another output of the same compound and the mean pairwise Tanimoto of the top 5 (`tmp/eos8zvb/nested_stages.py`).
- **Slurm tasks run the clone's code when each task starts:** cancel the array before editing, pulling or checking out the model clone, and archive its partial `results/<id>/{csv,logs,out}` (to `tmp/`) before re-submitting; otherwise the array mixes old and new code.
- **Pending tasks can wait a long time:** the lab CPU nodes are shared and were nearly full; a short array (eos8zvb tasks take about 20 s each on a node) may sit pending behind longer ones.

## Open at the time of writing

- 18 models not yet run on the benchmark; all their issue rows say "pending".
- 5 models still have Output Consistency = Fixed; eos694w has a red image test; eos935d's description says up to 10 metabolites while its Output Dimension is 15; 14 models have unreleased output-affecting changes.
- eos2401: the GPU and serial decode/filter change is a proposal only; decision pending.
- The chemical-space map and grid are drawn for eos2401 and eos8zvb (reference projection cached in `cache/eos1klk/`; older grids exist for eos4qda and eos9taz); draw the other models' with `07_plot_examples.py --model <id>`. When a new fragment model's outputs all grow from one fragment of the input, add it to `SHARED_FRAGMENT_MODELS` in script 07. Models that were not benchmarked use their shipped example output.
- Scripts 01-07, `gen-models-master.csv` and this file are committed on `chemsampler-update`; check `git status` for what has been pushed.
