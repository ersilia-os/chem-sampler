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

ChemSampler relies on the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia) to fetch and run generative models. To also run the test suite, install with the `dev` extra instead: `pip install -e ".[dev]"`.

## Usage

Generators and a similarity model propose molecules, and annotators score them.
`hill_climb` alternates the two, promoting each round's best candidate to seed the
next round. An annotator is either **primary** (optimized) or **secondary** (a
guard-rail: a candidate that does not satisfy it never wins, in any mode). How the
primaries combine is chosen with `mode`, which has no default. A lighter
`annotate(smiles, annotators)` scores one molecule with no generation or rounds,
which is handy for checking a seed's baseline:

- `mode="sequential"`: primaries are optimized one at a time, in list order.
  Each finished stage's own cutoff becomes a hard floor for every later stage,
  so a later stage can trade within an already-cleared objective's accepted
  range but never push it back below its declared goal. (If a stage's round
  budget runs out before reaching its own cutoff, its actual achieved value
  becomes the floor instead.)
- `mode="joint"`: all primaries optimized together; a candidate's score is its
  count of cutoffs satisfied, with ties broken arbitrarily.
- `mode="weighted"`: all primaries optimized together as one scalarized score,
  `sum(weight * value)` per annotator (negated first for `direction="lower"`
  annotators, so maximizing the sum always pushes every annotator the right
  way). `weight` defaults to 1.0. No cross-annotator normalization is applied,
  so pick weight magnitudes relative to each annotator's own scale (e.g. a 0-1
  QED score needs a much larger weight than a score in the hundreds to have
  comparable pull).
- `mode="incremental"`: like `sequential`, but each primary is advanced a
  little at a time instead of all at once, cycling through the primaries
  `n_steps` times. Each visit targets an interim value that interpolates
  linearly between the seed's own score and that annotator's cutoff (e.g. a
  seed at MAIP=10 with a cutoff of 300 and `n_steps=10` asks for MAIP>=39 on
  the first visit, MAIP>=68 on the second, ..., MAIP>=300 on the tenth). Once
  reached, an interim target becomes a floor for later steps, same as a
  finished `sequential` stage. Requires `seed_smiles`, `n_steps`, and a
  finite cutoff for every primary. Slower than `sequential`, useful for
  difficult cases where committing one annotator all the way to its cutoff
  before touching the next risks an unrecoverable overshoot.

```python
from chemsampler.config import load_inputs
from chemsampler.optimize import hill_climb

# reserpine, PubChem CID 5770
seed = "CO[C@H]1[C@@H](C[C@@H]2CN3CCC4=C([C@H]3C[C@@H]2[C@@H]1C(=O)OC)NC5=C4C=CC(=C5)OC)OC(=O)C6=CC(=C(C(=C6)OC)OC)OC"

inputs = load_inputs(
    path_to_models="examples/inputs/csv",
    path_to_model_repositories="/path/to/model/repositories",
)
summary, candidates_by_round = hill_climb(
    inputs.generators,
    inputs.annotators,
    similarity=inputs.similarity,
    mode="sequential",
    seed_smiles=seed,
    tanimoto_cutoff=0.4,
)
```

`summary` has one row per round (`round`, `smiles`, `score`, `is_new_best`,
`active_annotator_id` — which primary that round was driving, `None` in
`mode="joint"`). `candidates_by_round[n]` has one row per candidate considered
in round `n`, with its source model, every annotator's value,
`cutoffs_satisfied` (primaries only), `guardrails_satisfied` (only candidates
where it is true can win), `tanimoto_to_seed` (only if a seed was given), and
`tanimoto_to_original_seed` (only if `original_seed_smiles` was given).
`tanimoto_cutoff`/`tanimoto_direction` gate on `tanimoto_to_seed`:
`direction="higher"` (default) keeps candidates at least this similar to the
seed; `"lower"` pushes toward novelty instead.
`seed_smiles` is optional; without one, round 1 is unconditionally the new
best, but every Hub model needs a seed to run, so pass one in practice.

### Inputs

Annotators are required, and so is at least one of similarity and generators.
Each input comes from the first of these that has it: its own flag (a CSV file);
the config file's key of that name, or a CSV of that name in `--path-to-models`
(both is an error); a CSV of that name in the working directory, which is only
looked at when no folder is named. See [`examples/inputs/`](examples/inputs/) for
one model of each kind.

| Input | Flag (a CSV file) | CSV in `--path-to-models` | Config key |
|---|---|---|---|
| similarity | `--similarity` | `similarity.csv`, header `similarity`, one model (`eos8nyg` for now) | `"similarity": ["eos8nyg"]` |
| generators | `--generators` | `generators.csv`, header `generators` | `"generators": ["eos9taz", ...]` |
| annotators | `--annotators` | `annotators.csv`, header `source,model,cutoff,direction,priority` (+ `weight,column`) | `"annotators": [{...}]` |

An annotator row: `source` is `ersilia` (`model` is a Hub id) or `rdkit` (`model`
is `mw` or `qed`); `cutoff` is a number, or `min-max` with `direction` `range`
(`higher`/`lower` take a single number, and only a secondary can be a range);
`priority` is `primary` (optimized, in list order for `sequential` and
`incremental`) or `secondary` (a guard-rail, no weight); `weight` (primary,
default 1) and `column` (a multi-output Hub model) are optional.

```json
{
  "path_to_model_repositories": "/path/to/model/repositories",
  "similarity": ["eos8nyg"],
  "generators": ["eos9taz", "eos8vud"],
  "annotators": [
    {"source": "ersilia", "model": "eos4zfy", "cutoff": 0.0, "direction": "higher", "priority": "primary"},
    {"source": "rdkit", "model": "mw", "cutoff": "250-450", "direction": "range", "priority": "secondary"}
  ]
}
```

Paths in the config file are relative to the file; `--path-to-models` and
`--path-to-model-repositories` override its `path_to_models` and
`path_to_model_repositories`.

Every Hub model is checked against the catalog in
[`src/chemsampler/hub/catalog.csv`](src/chemsampler/hub/catalog.csv): generators
must be labelled Sampling / Generation, the similarity model Sampling /
Similarity search, annotators Annotation. The catalog lists only models the Hub
lists (the Ready ones), and a similarity model other than `eos8nyg` is refused.
The check is instant and works offline. `chemsampler setup` rebuilds the catalog
from the Hub (about 20 s, needs a network connection), tells you which models are
new or gone, and leaves the file untouched if the Hub cannot be reached. Commit
the file after running it.

A similarity-search model returns the seed's nearest neighbours in a reference
library: known compounds, not generated ones. Its neighbours of the current best
molecule join each round's candidates, and the `source` column shows
`eos8nyg`.

### Backends

Two backends are available for running models:

**`run_sh` (default)**: Shells out directly to a model's bundled `run.sh` script.
- **Requires:** Model folder at `<path-to-model-repositories>/{model_id}`
- **No:** persistent server, ports, or Ersilia initialization
- **Fast:** minimal startup overhead
- **Limitation:** only works for conda-packed models fetched locally

**`ersilia`**: Uses Ersilia's HTTP API (`ersilia serve`/`run`/`close`).
- **Requires:** Ersilia installed; model fetched via `ersilia fetch`
- **Works:** with any Ersilia model (Docker, conda, or other packaging)
- **Cost:** slower startup (server initialization) and process management overhead

Use `run_sh` by default. Use `ersilia` only for models not conda-packed or not
locally available (`--backend ersilia`; the repositories folder is then not needed).

## CLI

| Command | Description |
|---|---|
| `chemsampler run` | Generate and score candidate molecules across rounds, saving the result. |
| `chemsampler annotate` | Score a single molecule against a set of annotators (or one model), no generation. |
| `chemsampler setup` | Rebuild the model catalog from the Ersilia Model Hub and report what changed. |

```bash
chemsampler run \
  --path-to-models examples/inputs/csv \
  --path-to-model-repositories /path/to/model/repositories \
  --mode sequential \
  --seed-smiles "CO[C@H]1[C@@H](C[C@@H]2CN3CCC4=C([C@H]3C[C@@H]2[C@@H]1C(=O)OC)NC5=C4C=CC(=C5)OC)OC(=O)C6=CC(=C(C(=C6)OC)OC)OC" \
  --output-dir results/

chemsampler run --config config.json --mode joint --seed-smiles "CCO" --output-dir results/
chemsampler annotate --path-to-models examples/inputs/csv --smiles "CCO"
chemsampler annotate --source rdkit --model mw --smiles "CCO"
chemsampler setup
```

`--backend` picks between `"ersilia"` and `"run_sh"` for every Hub model in the
run — see [Backends](#backends) above. `annotate` needs only the annotators.
(The audit scripts in `examples/gen-models-audit/` use `--path-to-models` for the
folder of model clones, which is a different folder.)
Run `chemsampler run --help`/`chemsampler annotate --help` for the full
option lists.

## About the Ersilia Open Source Initiative

The [Ersilia Open Source Initiative](https://ersilia.io) is a tech-nonprofit organization fueling sustainable research in the Global South. Ersilia's main asset is the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia), an open-source repository of AI/ML models for antimicrobial drug discovery.

![Ersilia Logo](assets/Ersilia_Brand.png)
