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
chemsampler init -p my-run && cd my-run   # a project: input/ and runs/
chemsampler fetch eos9p57 eos8nyg eos4zfy rdkit # these models, live (writes input/generators.csv and similarity.csv), and the rdkit env
$EDITOR input/annotators.csv              # cutoff, direction and priority of each annotator
printf 'smiles\nCCO\n' > input/smiles.csv  # the seeds, one per row
chemsampler run                           # --mode only with several primaries
```

`init` only creates a new folder; without `-p` it makes a dummy `chemsampler-test`.

```
my-run/                       ~/.chemsampler/        the store ($CHEMSAMPLER_HOME)
  chemsampler.json              models/              one git clone per model
  input/                        envs/                one conda environment per model
  runs/                         logs/                what fetch ran, per model
                                hpc.json             cluster profiles, if any
```

Every location is set once, by `init`, and written to the project's
`chemsampler.json`; the other commands need only the project. They work anywhere
inside it, or on any project with `-p PATH`. Each run writes a dated folder under
`runs/` with an `inputs.json` of the models and the options it ran with, and one
folder per seed (`seed1/`, `seed2/`, ...: the seeds are climbed one after
another) with its `summary.csv` and one `round<n>.csv` per round. The paths are absolute: a copied project is refused until its
`chemsampler.json` names its own `input/` and `runs/`.

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

Every model is used at a **version**, given next to it in `input/*.csv`: a
release tag of its GitHub repository (e.g. `v1.0.0`), or `live`, the latest
commit of its main branch (which need not be a release). `chemsampler fetch ID [ID ...]`
(spaces or commas; `ID` is `live`, `ID@TAG` a given release) tells
generators, similarity models and annotators apart by itself, and needs at least
one annotator and one generator or similarity model. It sets
`input/generators.csv` and `input/similarity.csv` to the ids and their versions,
asking before it changes them (`-f`: without asking); `annotators.csv` stays
yours (its versions are the ones fetched), and fetch says where it disagrees with
the ids. Without ids, it lists everything `input/*.csv` names, the `rdkit`
environment included, and fetches it once you agree (`-f`: without asking).
Either way it is safe to repeat, and leaves every model at its version:

- A missing clone is cloned and checked out at its release; an existing one is
  switched to it (`live`: to its main branch, then pulled). A clone with local changes or local commits is never touched,
  and `run` refuses a clone that is not at the version its row names. Checkpoints
  come from S3 as they are now: they are not versioned.
- An environment is built when missing, and rebuilt (after switching, so against
  that release's recipe) when a pinned package differs from its recipe (`install.yml`, or the
  Dockerfile for older models), e.g. `torch: installed 2.8.0, recipe 2.9.1`.
  Recipe lines that pin nothing (an unpinned package, a GitHub install, a shell
  line) cannot be compared and count as matching; `fetch` says how many. Don't
  fetch while a run uses the same models.
- An environment left half-built by an interrupted run is removed and rebuilt.
- The RDKit annotators (`mw`, `qed`, `tanimoto`) run in a small `rdkit`
  environment, here and on a cluster: `chemsampler fetch rdkit` builds it (alone,
  or with model ids), and so does plain `fetch` when `annotators.csv` lists one;
  `run` refuses to start without it.

## Inputs

Annotators are required, and so is at least one of similarity and generators.

| Input | File in `input/` | What it holds |
|---|---|---|
| similarity | `similarity.csv` | `model,version`: one similarity-search model, `eos8nyg` (Ersilia Reference Library), `eos1bgt` (COCONUT) or `eos80al` (Enamine Hit Locator) |
| generators | `generators.csv` | `model,version`: generative models |
| annotators | `annotators.csv` | models or RDKit functions that score molecules, each with a cutoff, a direction and a priority |
| seeds | `smiles.csv` | the starting molecules, header `smiles`, one per row; each is climbed in turn (`run --seed-smiles` gives one instead) |

A **primary** annotator is optimized; a **guard-rail** is one every candidate must
satisfy. Any annotator can be either, including the RDKit ones: `mw`, `qed` and
`tanimoto`, the similarity to the seed of each climb (e.g. `rdkit,tanimoto,0.4,higher,guard-rail`
keeps candidates near it). Every model is checked against a bundled catalog of what
the Hub currently labels for each role, which works offline. The exact columns, the
project file and the catalog are in [docs/inputs.md](docs/inputs.md). With one
primary a run is sequential; with several, `--mode` picks one of the four ways to
combine them (`sequential`, `joint`, `weighted`, `incremental`), described with the
Python API in [docs/modes.md](docs/modes.md).

## Running on a cluster

`chemsampler run --hpc NAME` sends every model of a round to a Slurm cluster as its
own job: the similarity and generative models side by side, then, once their
candidates are put together, one job per annotator (`mw`, `qed` and `tanimoto`
included). The terminal waits, checking every 10 s (every 60 s once a step passes
5 minutes), and says how many jobs are pending, running and done; each round is
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
| `chemsampler fetch` | Clone the models at their version, download their checkpoints and build their conda environments, in the store (asks first); `fetch rdkit` builds the RDKit annotators' one. |
| `chemsampler run` | Generate and score candidate molecules across rounds, saving the result. |
| `chemsampler resume` | Continue a stopped run from the rounds it finished. |
| `chemsampler annotate` | Score a single molecule against a set of annotators (or one model), no generation. |
| `chemsampler catalog` | List the models on this machine; `--hub` lists what the Hub offers, `--update` refreshes the bundled catalog. |
| `chemsampler store` | Show where models, environments and cluster profiles are kept, and how to move them. |

```bash
chemsampler fetch eos9p57,eos8nyg rdkit -f      # set the project's samplers (live) and fetch them
chemsampler fetch eos9p57@v1.0.0 eos8nyg rdkit   # a given release of eos9p57
chemsampler fetch                                # everything input/*.csv lists, after asking
chemsampler run -p my-run                        # from outside the project
chemsampler run --mode weighted                  # several primaries: say how to combine them
chemsampler run --seed-smiles "CCO"              # one seed instead of input/smiles.csv
chemsampler run --hpc mycluster                  # models as Slurm jobs
chemsampler run --hpc mycluster --detach         # the run is a job too
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
