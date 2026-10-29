# Modes and the Python API

Generators and a similarity model propose molecules, and annotators score them.
`hill_climb` alternates the two, promoting each round's best candidate to seed the
next round. An annotator is either **primary** (optimized) or **secondary** (a
guard-rail: a candidate that does not satisfy it never wins, in any mode). How the
primaries combine is chosen with `mode`, which has no default. A lighter
`annotate(smiles, annotators)` scores one molecule with no generation or rounds,
which is handy for checking a seed's baseline.

## Modes

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

## From Python

```python
from chemsampler.config import load_inputs
from chemsampler.optimize import hill_climb

# reserpine, PubChem CID 5770
seed = "CO[C@H]1[C@@H](C[C@@H]2CN3CCC4=C([C@H]3C[C@@H]2[C@@H]1C(=O)OC)NC5=C4C=CC(=C5)OC)OC(=O)C6=CC(=C(C(=C6)OC)OC)OC"

# run from inside a project: its inputs and the store's models are found for you
inputs = load_inputs()
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

`hill_climb`'s generator may be `None` when a similarity model is given. Models
the code does not know are rejected: generators must be labelled Sampling /
Generation in the catalog, the similarity model Sampling / Similarity search, and
annotators Annotation.
