# Sampling chemical space around a molecule

ChemSampler generates and ranks new candidate molecules around a seed compound, using generative models from the [Ersilia Model Hub](https://ersilia.io/model-hub). It is developed and maintained by the [Ersilia Open Source Initiative](https://ersilia.io).

> A legacy version of ChemSampler is available in the archived [`chem-sampler-legacy`](https://github.com/ersilia-os/chem-sampler-legacy) repository.

## Status

:construction: Early development. The pipeline below works, but the API is unstable and there is no CLI yet.

## Installation

```bash
conda create -n chemsampler python=3.10
conda activate chemsampler
pip install git+https://github.com/ersilia-os/chem-sampler.git
```

ChemSampler relies on the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia) to fetch and run generative models.

## Usage

A generator proposes molecules and a set of annotators scores them. `hill_climb`
alternates the two, promoting each round's best candidate to seed the next round.
Every annotator is uniform — a `cutoff` and a `direction` ("higher" or "lower" is
better), nothing distinguishes one from another. How they combine is chosen with
`mode`, which has no default:

- `mode="sequential"`: annotators are optimized one at a time, in list order.
  Each finished stage's achieved value becomes a hard floor for every later
  stage, so a later stage can never trade away an earlier gain.
- `mode="joint"`: all annotators optimized together; a candidate's score is its
  count of cutoffs satisfied, with ties broken arbitrarily.

```python
from chemsampler.models.annotator import QEDAnnotator
from chemsampler.models.generator import HubGenerator
from chemsampler.models.hub_annotator import HubAnnotator
from chemsampler.models.spec import AnnotatorSpec
from chemsampler.optimize import hill_climb

# reserpine, PubChem CID 5770
seed = "CO[C@H]1[C@@H](C[C@@H]2CN3CCC4=C([C@H]3C[C@@H]2[C@@H]1C(=O)OC)NC5=C4C=CC(=C5)OC)OC(=O)C6=CC(=C(C(=C6)OC)OC)OC"

summary, candidates_by_round = hill_climb(
    generator=HubGenerator("eos9taz"),
    annotators=[
        AnnotatorSpec("eos4zfy", HubAnnotator("eos4zfy"), cutoff=0.0, direction="higher"),  # MAIP
        AnnotatorSpec("qed", QEDAnnotator(), cutoff=0.3, direction="higher"),  # drug-likeness floor
    ],
    mode="sequential",
    seed_smiles=seed,
    tanimoto_cutoff=0.4,
)
```

`summary` has one row per round (`round`, `smiles`, `score`, `is_new_best`,
`active_annotator_id` — which annotator that round was driving, `None` in
`mode="joint"`). `candidates_by_round[n]` has one row per candidate considered
in round `n`, with its source generator, every annotator's value,
`cutoffs_satisfied`, and `tanimoto_to_seed` (only if a seed was given).
`seed_smiles` is optional; without one, round 1 is unconditionally the new
best, and only seed-agnostic generators (like `ChemblSampler`) can contribute
to it.

`HubGenerator` accepts any Hub model that returns generated molecules. Several
can be pooled behind the same interface, so a pool drops in wherever a single
generator is expected:

```python
from chemsampler.models.generator import GeneratorPool

pool = GeneratorPool([HubGenerator("eos9taz"), HubGenerator("eos6ost")])
```

`ChemblSampler` is a null baseline: it draws molecules at random from a filtered
ChEMBL set (single-component, 200-450 Da) and **ignores the seed**. It answers
"does a generator beat a random draw from known chemistry?". It is opt-in only —
never included by default:

```python
from chemsampler.models.chembl import ChemblSampler

baseline = ChemblSampler(n=1000, random_state=42)
```

The reference set is fetched with [`eosvc`](https://github.com/ersilia-os/eosvc) on
first use, or rebuilt from scratch with `python -m chemsampler.data.chembl`.

### Backends

`HubGenerator` and `HubAnnotator` accept `backend="run_sh"` as an alternative to
the default `backend="ersilia"`. The default serves each model over HTTP
(`ersilia serve`/`run`/`close`); `run_sh` instead shells out directly to the
model's bundled `run.sh`, with no persistent server — it avoids the
orphaned-process and port-contention issues the default backend can hit under
sustained use, but only works for a model that's already fetched locally and
conda-packed:

```python
HubGenerator("eos9taz", backend="run_sh")
```

### Config files

Generators and annotators can also be loaded from CSVs instead of built by hand:

```python
from chemsampler.config import load_annotators, load_generators

generators = load_generators()  # defaults to the 3 validated Hub generators
annotators = load_annotators("my_annotators.csv")
```

`load_generators(path=None)` reads a `generator_id` column; each id is either an
Ersilia model id, or the literal `"chembl"` to opt into `ChemblSampler` — the
shipped default never includes it, so it must be listed explicitly in your own
CSV. `load_annotators(path)` reads `annotator_id, cutoff, direction, column`
(`cutoff` and `direction` are required for every row; row order is the
priority order used by `mode="sequential"`; `column` is optional and picks an
output column for a multi-output Hub model).

See [`examples/`](examples/) for a runnable script.

## About the Ersilia Open Source Initiative

The [Ersilia Open Source Initiative](https://ersilia.io) is a tech-nonprofit organization fueling sustainable research in the Global South. Ersilia's main asset is the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia), an open-source repository of AI/ML models for antimicrobial drug discovery.

![Ersilia Logo](assets/Ersilia_Brand.png)
