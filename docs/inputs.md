# Inputs

Annotators are required, and so is at least one of similarity and generators.
Each input comes from the first of these that has it:

1. its own flag, a CSV file (`--similarity`, `--generators`, `--annotators`);
2. the config file's key of that name, or a CSV of that name in the folder given by
   `--path-to-models` (both is an error);
3. a CSV of that name in the working directory, which is only looked at when no
   folder is named.

Inside a project, its `chemsampler.json` is the config file and `input/` is the
folder. See [`examples/inputs/`](../examples/inputs/) for one model of each kind.

| Input | Flag (a CSV file) | CSV in `input/` | Config key |
|---|---|---|---|
| similarity | `--similarity` | `similarity.csv`, header `similarity`, one model (`eos8nyg` for now) | `"similarity": ["eos8nyg"]` |
| generators | `--generators` | `generators.csv`, header `generators` | `"generators": ["eos9taz", ...]` |
| annotators | `--annotators` | `annotators.csv`, header `source,model,cutoff,direction,priority` (+ `weight,column`) | `"annotators": [{...}]` |

## Annotators

`source` is `ersilia` (`model` is a Hub id) or `rdkit` (`model` is `mw` or `qed`).
`cutoff` is a number, or `min-max` with `direction` `range` (`higher`/`lower` take a
single number, and only a secondary can be a range). `priority` is `primary`
(optimized, in list order for `sequential` and `incremental`) or `secondary` (a
guard-rail, no weight). `weight` (primary, default 1) and `column` (a multi-output
Hub model) are optional.

## A config file

```json
{
  "similarity": ["eos8nyg"],
  "generators": ["eos9taz", "eos8vud"],
  "annotators": [
    {"source": "ersilia", "model": "eos4zfy", "cutoff": 0.0, "direction": "higher", "priority": "primary"},
    {"source": "rdkit", "model": "mw", "cutoff": "250-450", "direction": "range", "priority": "secondary"}
  ]
}
```

Besides the three lists it may name `path_to_models`, `path_to_rounds`,
`path_to_model_repositories` and `path_to_envs`. Paths in it are relative to the
file; the flags override it. A project that names its own
`path_to_model_repositories` and `path_to_envs` keeps its models apart from the
store.

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
