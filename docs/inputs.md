# Inputs

Annotators are required, and so is at least one of similarity and generators.
They are three CSV files in the project's `input/` folder, with an optional fourth
for the seed; `chemsampler init`
writes a template of each (the same as [`examples/inputs/csv/`](../examples/inputs/csv/)).
`chemsampler fetch ID [ID ...]` sets `similarity.csv` and `generators.csv` to the
ids it is given (a role with none has its file removed); `annotators.csv` is
always yours to write.

| Input | CSV in `input/` |
|---|---|
| similarity | `similarity.csv`, header `similarity`, one model: `eos8nyg`, `eos1bgt` or `eos80al` |
| generators | `generators.csv`, header `generators` |
| annotators | `annotators.csv`, header `source,model,cutoff,direction,priority` (+ `weight,column`) |
| seed | `smiles.csv`, header `smiles`, one molecule (optional) |

`run` takes its seed from `smiles.csv`; `--seed-smiles` overrides it, with a note.
For now the file holds exactly one molecule
([several seeds](https://github.com/ersilia-os/chem-sampler/issues/8) are still to
be defined); without the file a run has no seed. `init` writes no `smiles.csv`.

## Annotators

`source` is `ersilia` (`model` is a Hub id) or `rdkit` (`model` is `mw`, `qed` or
`tanimoto`). `cutoff` is a number, or `min-max` with `direction` `range`
(`higher`/`lower` take a single number, and only a guard-rail can be a range).
`priority` is `primary` (optimized, in list order for `sequential` and
`incremental`) or `guard-rail` (a constraint every winner must satisfy, no weight).

`tanimoto` is the similarity of each candidate to the run's seed (Morgan
fingerprints, radius 2, 2048 bits; the seed itself scores 1.0), so it needs a seed.
As a guard-rail it keeps candidates near the seed (`0.4,higher`), away from it
(`0.6,lower`) or both (`0.4-0.9,range`); as a primary, `lower` pushes the search
away from the seed. The run computes it itself: on a cluster it is no job and
needs no `rdkit` environment. `weight` (primary, default 1) and `column` (a multi-output
Hub model) are optional.

## The project file

`chemsampler init` writes every location into the project's `chemsampler.json`,
as absolute paths:

```json
{
  "path_to_models": "/home/you/runs/my-run/input",
  "path_to_rounds": "/home/you/runs/my-run/rounds",
  "path_to_model_repositories": "/home/you/.chemsampler/models",
  "path_to_envs": "/home/you/.chemsampler/envs"
}
```

`path_to_models` is the folder with the three CSVs, `path_to_rounds` where runs
write their results, and the last two where `fetch` puts the model clones and
their conda environments: this machine's store unless `init --models DIR
--envs DIR` named others. To change a location later, edit the file (a relative
path there is relative to the project folder); never move model or environment
folders after `fetch`.

A copied or moved project still names the old folders. Commands refuse a project
whose input or rounds folder belongs to another project, so a copy never reads the
original's inputs or writes into its results: fix the two paths in the copy's
`chemsampler.json`.

## The catalog

Every Hub model is checked against
[`src/chemsampler/hub/catalog.csv`](../src/chemsampler/hub/catalog.csv): generators
must be labelled Sampling / Generation, the similarity model Sampling / Similarity
search, annotators Annotation. The catalog lists only models the Hub lists (the
Ready ones). Of the Hub's similarity models only three are accepted, the ones that
run the same search on a fixed local library: `eos8nyg` (Ersilia Reference
Library), `eos1bgt` (COCONUT natural products) and `eos80al` (Enamine Hit
Locator); the others call outside services. The check is instant and works
offline.

`chemsampler catalog --update` rebuilds the catalog from the Hub (about 20 s, needs
a network connection), tells you which models are new or gone, and leaves the file
untouched if the Hub cannot be reached. Commit the file after running it.

A similarity-search model returns the seed's nearest neighbours in a reference
library: known compounds, not generated ones. Its neighbours of the current best
molecule join each round's candidates, and the `source` column shows its id.
