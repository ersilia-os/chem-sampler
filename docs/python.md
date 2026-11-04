# How a run works, and the Python API

Generators and similarity models propose molecules, and annotators score them.
`hill_climb` alternates the two, promoting each round's best candidate to seed the
next round. Exactly one annotator is **primary** (optimized); every other one is
a **guard-rail** (a candidate that does not satisfy it never wins). A lighter
`annotate(smiles, annotators, seed_smiles)` scores one molecule with no generation
or rounds, which is handy for checking a seed's baseline.

## How a round picks its winner

1. Every candidate the models propose is scored by every annotator.
2. Candidates that fail any guard-rail are set aside; they can never win.
3. Among the rest, the best primary value wins: the highest for `higher`, the
   lowest for `lower`. A candidate the primary could not score never wins.
4. The winner becomes the new best, and the next round's starting molecule, if
   it beats the current best by more than `tolerance`.

The run stops when a round does not improve, when no candidate passes the
guard-rails (the best stays the previous round's, or the seed), when the
primary's cutoff is met, or after `n_rounds`. A primary whose cutoff is `nan`
has none: it is optimized until a round does not improve or the rounds run out.
An infinite cutoff never stops a run early either.

## From Python

```python
from chemsampler.config import load_inputs
from chemsampler.optimize import hill_climb

# run from inside a project: its inputs, its seeds (input/smiles.csv) and the
# store's models are found for you; this climbs from the first seed
inputs = load_inputs()
summary, candidates_by_round = hill_climb(
    inputs.generators,
    inputs.annotators,
    similarity=inputs.similarity,
    seed_smiles=inputs.seeds[0],
)
```

`summary` has one row per round (`round`, `smiles`, `score`, `is_new_best`).
`candidates_by_round[n]` has one row per candidate considered in round `n`, best
first, with its source model, every annotator's value, `cutoff_satisfied` (the
primary's; empty when it has no cutoff), `guardrails_satisfied` (only candidates
where it is true can win) and `tanimoto_to_seed` (only if a seed was given).
To limit similarity to the seed, add a Tanimoto annotator, e.g.
`AnnotatorSpec("tanimoto", TanimotoAnnotator(), 0.4, "higher", priority="guard-rail")`
(from `chemsampler.models.annotator`); `hill_climb` compares it with
`seed_smiles`. `"lower"` pushes toward novelty instead, and as the primary it is
optimized like any annotator. RDKit annotators built by hand, as here, are
computed with chemsampler's own RDKit; those `load_inputs` builds run in the
project's `rdkit` environment (`chemsampler fetch rdkit`), with the same values.
`seed_smiles` is optional (but required with a Tanimoto annotator); without one, round 1 is unconditionally the new
best, but every Hub model needs a seed to run, so pass one in practice.

`hill_climb`'s generator may be `None` when a similarity model is given, and
`similarity` takes one model or a list. Models
the code does not know are rejected: generators must be labelled Sampling /
Generation in the catalog, the similarity models Sampling / Similarity search, and
annotators Annotation.
