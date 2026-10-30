# Sampling chemical space around a molecule

ChemSampler generates and ranks new candidate molecules around a seed compound, using generative models from the [Ersilia Model Hub](https://ersilia.io/model-hub). It is developed and maintained by the [Ersilia Open Source Initiative](https://ersilia.io).

> A legacy version of ChemSampler is available in the archived [`chem-sampler-legacy`](https://github.com/ersilia-os/chem-sampler-legacy) repository.

## Status

:construction: Early development. The pipeline below works, but the API is unstable.

## Installation

Clone the repository first, then create the environment from inside the clone (the editable install points at it):

```bash
git clone https://github.com/ersilia-os/chem-sampler.git
cd chem-sampler
conda create -n chemsampler python=3.10
conda activate chemsampler
pip install -e .
```

To also run the test suite, install with the `dev` extra: `pip install -e ".[dev]"`.

## Getting started

An analysis is a **project**: a folder with its inputs and results. The models it
uses are expensive to fetch, so they live in a **store** that every project on the
machine shares.

```bash
chemsampler init -p my-run && cd my-run   # a project: input/ and rounds/
$EDITOR input/*.csv                       # which models to use (templates are written for you)
chemsampler fetch                         # clone them, download checkpoints, build their conda envs
chemsampler run --mode sequential --seed-smiles "CCO"
```

`init` only creates a new folder; without `-p` it makes a dummy `chemsampler-test`.

```
my-run/                       ~/.chemsampler/        the store
  chemsampler.json              models/              one git clone per model
  input/                        envs/                one conda environment per model
  rounds/                       logs/                what fetch ran, per model
```

Commands work anywhere inside a project and take their folders from its
`chemsampler.json`, so `run` needs no paths and no `--output-dir`: each run writes
a dated folder under `rounds/` with its `summary.csv`, one `round<n>.csv` per round
and an `inputs.json` of what was resolved. Switching project is `cd`.

**Where models live.** A second project reuses what the first one fetched, so
`fetch` there is instant. `chemsampler catalog` shows what the store holds and the
state of each environment. To keep the store somewhere else (on a cluster, a
shared filesystem), write `~/.chemsampler/config.json`:

```json
{"store": "/shared/chemsampler"}
```

or move just one folder with `path_to_model_repositories` / `path_to_envs`, or name
them with `--path-to-model-repositories` / `--path-to-envs`, or in a project's
`chemsampler.json` (which then keeps its models apart). Two cautions: **never
move or rename a store** after `fetch`, because conda environments record their own
path; and prefer a local disk, since model start-up was several times slower on a
shared filesystem in our tests.

`fetch` is safe to repeat. It clones only what is missing and reports the rest:

- An existing clone is only inspected (`up to date`, `3 behind`, `local changes`).
  `fetch --update` pulls the ones that are behind, fast-forward only, and never
  touches a clone with local changes or local commits.
- An existing environment is compared with its model's recipe (`install.yml`, or
  the Dockerfile for older models) and differences are reported, for example
  `torch: installed 2.8.0, recipe 2.9.1`. `fetch --rebuild` rebuilds only the
  environments that differ. Recipe lines that pin nothing (an unpinned package, a
  GitHub install, a shell line) cannot be compared, and `fetch` says how many it
  skipped.
- An environment left half-built by an interrupted run is removed and rebuilt.

## Inputs

Annotators are required, and so is at least one of similarity and generators.

| Input | File in `input/` | What it holds |
|---|---|---|
| similarity | `similarity.csv` | one similarity-search model (`eos8nyg` for now) |
| generators | `generators.csv` | generative models |
| annotators | `annotators.csv` | models or RDKit functions that score molecules, each with a cutoff, a direction and a priority |

A **primary** annotator is optimized; a **secondary** one is a guard-rail that every
candidate must satisfy. Every model is checked against a bundled catalog of what
the Hub currently labels for each role, which works offline. Flags, a config file,
the exact columns and the catalog are in [docs/inputs.md](docs/inputs.md); the four
ways to combine primaries (`sequential`, `joint`, `weighted`, `incremental`) and the
Python API are in [docs/modes.md](docs/modes.md).

## Backends

`run_sh` (default) runs a model's `run.sh` directly: fast, no server, but it needs
the model fetched locally (what `chemsampler fetch` does, or `ersilia fetch`).
`ersilia` (`--backend ersilia`) goes through `ersilia serve`/`run`: slower to start,
but it works with any Ersilia model and needs no store.

## CLI

| Command | Description |
|---|---|
| `chemsampler init` | Create a project for one analysis's inputs and results. |
| `chemsampler fetch` | Clone the models, download their checkpoints and build their conda environments, in the store. |
| `chemsampler run` | Generate and score candidate molecules across rounds, saving the result. |
| `chemsampler annotate` | Score a single molecule against a set of annotators (or one model), no generation. |
| `chemsampler catalog` | List the models on this machine; `--hub` lists what the Hub offers, `--update` refreshes the bundled catalog. |

```bash
chemsampler fetch --update --rebuild   # pull clones that are behind, rebuild envs that differ
chemsampler annotate --smiles "CCO"
chemsampler annotate --source rdkit --model mw --smiles "CCO"
chemsampler catalog                    # * marks what the current project needs
chemsampler catalog --hub --role generators

# outside a project, say where everything is
chemsampler run --config config.json --mode joint --seed-smiles "CCO" --output-dir results/
```

(The audit scripts in `examples/gen-models-audit/` use `--path-to-models` for their
own folder of model clones, which is a different folder.) Run
`chemsampler <command> --help` for the full option lists.

## About the Ersilia Open Source Initiative

The [Ersilia Open Source Initiative](https://ersilia.io) is a tech-nonprofit organization fueling sustainable research in the Global South. Ersilia's main asset is the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia), an open-source repository of AI/ML models for antimicrobial drug discovery.

![Ersilia Logo](assets/Ersilia_Brand.png)
