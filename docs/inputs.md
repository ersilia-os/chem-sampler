# Inputs

Annotators are required, and so is at least one of similarity and generators.
They are three CSV files in the project's `input/` folder; `chemsampler init`
writes a template of each (the same as [`examples/inputs/csv/`](../examples/inputs/csv/)).

| Input | CSV in `input/` |
|---|---|
| similarity | `similarity.csv`, header `similarity`, one model (`eos8nyg` for now) |
| generators | `generators.csv`, header `generators` |
| annotators | `annotators.csv`, header `source,model,cutoff,direction,priority` (+ `weight,column`) |

## Annotators

`source` is `ersilia` (`model` is a Hub id) or `rdkit` (`model` is `mw` or `qed`).
`cutoff` is a number, or `min-max` with `direction` `range` (`higher`/`lower` take a
single number, and only a secondary can be a range). `priority` is `primary`
(optimized, in list order for `sequential` and `incremental`) or `secondary` (a
guard-rail, no weight). `weight` (primary, default 1) and `column` (a multi-output
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
Ready ones), and a similarity model other than `eos8nyg` is refused. The check is
instant and works offline.

`chemsampler catalog --update` rebuilds the catalog from the Hub (about 20 s, needs
a network connection), tells you which models are new or gone, and leaves the file
untouched if the Hub cannot be reached. Commit the file after running it.

A similarity-search model returns the seed's nearest neighbours in a reference
library: known compounds, not generated ones. Its neighbours of the current best
molecule join each round's candidates, and the `source` column shows `eos8nyg`.
