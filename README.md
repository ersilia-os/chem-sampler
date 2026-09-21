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

A generator proposes molecules and an annotator scores them. `hill_climb` alternates the two, promoting the best-scoring candidate to seed the next round, and stops once a round fails to improve.

```python
from chemsampler.models.annotator import QEDAnnotator
from chemsampler.models.generator import HubGenerator
from chemsampler.optimize import hill_climb

# reserpine, PubChem CID 5770
seed = "CO[C@H]1[C@@H](C[C@@H]2CN3CCC4=C([C@H]3C[C@@H]2[C@@H]1C(=O)OC)NC5=C4C=CC(=C5)OC)OC(=O)C6=CC(=C(C(=C6)OC)OC)OC"

df = hill_climb(
    seed_smiles=seed,
    generator=HubGenerator("eos9taz"),
    annotator=QEDAnnotator(),
)
```

The result has one row per round, with columns `round`, `smiles`, `score` and `is_new_best`.

`HubGenerator` accepts any Hub model that returns generated molecules. Several can be pooled behind the same interface, so a pool drops in wherever a single generator is expected:

```python
from chemsampler.models.generator import GeneratorPool, VALIDATED_GENERATORS

pool = GeneratorPool([HubGenerator(model_id) for model_id in VALIDATED_GENERATORS])
```

`ChemblSampler` is a null baseline: it draws molecules at random from a filtered
ChEMBL set (single-component, 200-450 Da) and **ignores the seed**. It answers
"does a generator beat a random draw from known chemistry?". It is opt-in, and
deliberately absent from `VALIDATED_GENERATORS`:

```python
from chemsampler.models.chembl import ChemblSampler

baseline = ChemblSampler(n=1000, random_state=42)
```

The reference set is fetched with [`eosvc`](https://github.com/ersilia-os/eosvc) on
first use, or rebuilt from scratch with `python -m chemsampler.data.chembl`.

See [`examples/`](examples/) for a runnable script.

## About the Ersilia Open Source Initiative

The [Ersilia Open Source Initiative](https://ersilia.io) is a tech-nonprofit organization fueling sustainable research in the Global South. Ersilia's main asset is the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia), an open-source repository of AI/ML models for antimicrobial drug discovery.

![Ersilia Logo](assets/Ersilia_Brand.png)
