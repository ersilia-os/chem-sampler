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

An analysis is a **project**: a folder with its inputs, its results and a note of
where its models are. The models are expensive to fetch, so by default they live in
a **store** that every project on the machine shares.

```bash
chemsampler init -p my-run && cd my-run   # a project: input/ and rounds/
$EDITOR input/*.csv                       # which models to use (templates are written for you)
chemsampler fetch                         # clone them, download checkpoints, build their conda envs
chemsampler run --mode sequential --seed-smiles "CCO"
```

`init` only creates a new folder; without `-p` it makes a dummy `chemsampler-test`.

```
my-run/                       ~/.chemsampler/        the store ($CHEMSAMPLER_HOME)
  chemsampler.json              models/              one git clone per model
  input/                        envs/                one conda environment per model
  rounds/                       logs/                what fetch ran, per model
                                hpc.json             cluster profiles, if any
```

Every location is set once, by `init`, and written to the project's
`chemsampler.json`; the other commands need only the project. They work anywhere
inside it, or on any project with `-p PATH`. Each run writes a dated folder under
`rounds/` with its `summary.csv`, one `round<n>.csv` per round and an `inputs.json`
of the models and the options it ran with. The paths are absolute: a copied project is refused until its
`chemsampler.json` names its own `input/` and `rounds/`.

Rounds are saved as they finish, with a copy of the inputs. A run stopped in any
way (Ctrl-C, a failed model, a closed terminal, a cancelled job) continues with
`chemsampler resume`: the latest run, or `resume RUN` for another. Its finished
rounds are read back, not run again, and the round it was in is run again whole.

**Where models live.** A second project reuses what the first one fetched, so
`fetch` there is instant. `chemsampler catalog` shows what the store holds and the
state of each environment. To put a project's models elsewhere (on a cluster, a
shared filesystem), give the folders to `init`:

```bash
chemsampler init -p my-run --models /shared/chemsampler/models --envs /local/envs
```

To move the store itself, and so the default of every new project, set
`CHEMSAMPLER_HOME` (e.g. `export CHEMSAMPLER_HOME=/shared/.chemsampler` in
`~/.bashrc`); everything chemsampler keeps is in that one folder. To change a
project later, edit its `chemsampler.json`. Two cautions: **never move or rename these folders** after
`fetch`, because conda environments record their own path; and prefer a local disk
for environments, since model start-up was several times slower on a shared
filesystem in our tests (a cluster run needs them on the shared one, see below).

`fetch` is safe to repeat. It clones only what is missing and reports the rest:

- An existing clone is only inspected (`up to date`, `3 behind`, `local changes`).
  `fetch --update` pulls the ones that are behind, fast-forward only, and never
  touches a clone with local changes or local commits.
- An existing environment is compared with its model's recipe (`install.yml`, or
  the Dockerfile for older models) and differences are reported, for example
  `torch: installed 2.8.0, recipe 2.9.1`. `fetch --update` also rebuilds the
  environments that differ, after pulling, so against the new recipe. Recipe lines that pin nothing (an unpinned package, a
  GitHub install, a shell line) cannot be compared, and `fetch` says how many it
  skipped.
- An environment left half-built by an interrupted run is removed and rebuilt.
- With an RDKit annotator in the inputs, `fetch` also builds a small `rdkit`
  environment, which a cluster run scores it with.

## Inputs

Annotators are required, and so is at least one of similarity and generators.

| Input | File in `input/` | What it holds |
|---|---|---|
| similarity | `similarity.csv` | one similarity-search model: `eos8nyg` (Ersilia Reference Library), `eos1bgt` (COCONUT) or `eos80al` (Enamine Hit Locator) |
| generators | `generators.csv` | generative models |
| annotators | `annotators.csv` | models or RDKit functions that score molecules, each with a cutoff, a direction and a priority |

A **primary** annotator is optimized; a **secondary** one is a guard-rail that every
candidate must satisfy. Every model is checked against a bundled catalog of what
the Hub currently labels for each role, which works offline. The exact columns, the
project file and the catalog are in [docs/inputs.md](docs/inputs.md); the four
ways to combine primaries (`sequential`, `joint`, `weighted`, `incremental`) and the
Python API are in [docs/modes.md](docs/modes.md).

## Running on a cluster

`chemsampler run --hpc NAME` sends every model of a round to a Slurm cluster as its
own job: the similarity and generative models side by side, then, once their
candidates are put together, one job per annotator (RDKit ones included). The
terminal waits and says how many jobs are pending, running and done; each round is
ranked as a local run would rank it. Ctrl-C cancels the run's jobs and keeps the
rounds already finished.

The cluster is described once, in `hpc.json` in the store, never in a project:

```json
{"mycluster": {"ssh": "login.example.org", "partition": "cpu",
               "nodes": ["node01", "node02"], "cpus": 4, "mem": "8G", "time": "01:00:00"}}
```

The project and the store must be on a filesystem the nodes share, at the same
paths (set `CHEMSAMPLER_HOME` there); start long runs inside `tmux`. Or add
`--detach`: the run itself becomes a small job (`env` and `driver` in the
profile), writes what it would print to `run.log` in its folder, and resumes by
itself if the cluster preempts it. Every key, and what a run
leaves in its folder, are in [docs/hpc.md](docs/hpc.md).

## CLI

| Command | Description |
|---|---|
| `chemsampler init` | Create a project for one analysis's inputs and results. |
| `chemsampler fetch` | Clone the models, download their checkpoints and build their conda environments, in the store. |
| `chemsampler run` | Generate and score candidate molecules across rounds, saving the result. |
| `chemsampler resume` | Continue a stopped run from the rounds it finished. |
| `chemsampler annotate` | Score a single molecule against a set of annotators (or one model), no generation. |
| `chemsampler catalog` | List the models on this machine; `--hub` lists what the Hub offers, `--update` refreshes the bundled catalog. |
| `chemsampler store` | Show where models, environments and cluster profiles are kept, and how to move them. |

```bash
chemsampler fetch --update                       # pull clones that are behind, rebuild envs that differ
chemsampler run -p my-run --mode joint --seed-smiles "CCO"   # from outside the project
chemsampler run --mode weighted --seed-smiles "CCO" --hpc mycluster   # models as Slurm jobs
chemsampler run --mode weighted --seed-smiles "CCO" --hpc mycluster --detach   # the run is a job too
chemsampler resume                               # continue the latest run where it stopped
chemsampler annotate --smiles "CCO"
chemsampler annotate --source rdkit --model mw --smiles "CCO"   # one model, no project needed
chemsampler catalog                              # * marks what the current project needs
chemsampler catalog --hub --role generators
```

(The audit scripts in `examples/gen-models-audit/` have their own `--path-to-models`
flag for their own folder of model clones, unrelated to these commands.) Run
`chemsampler <command> --help` for the full option lists.

## About the Ersilia Open Source Initiative

The [Ersilia Open Source Initiative](https://ersilia.io) is a tech-nonprofit organization fueling sustainable research in the Global South. Ersilia's main asset is the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia), an open-source repository of AI/ML models for antimicrobial drug discovery.

![Ersilia Logo](assets/Ersilia_Brand.png)
