# Inputs

Annotators are required, and so is at least one of similarity and generators.
They are three CSV files in the project's `input/` folder, and a fourth for the
seeds; `chemsampler init` writes a template of the first three (one model of
each kind, `live`). `chemsampler fetch ID [ID ...]` sets `similarity.csv` and
`generators.csv` to the ids it is given, `live` or at the release of `ID@TAG` (a
role with none has its file removed); `annotators.csv` is always yours to write.

| Input | CSV in `input/` |
|---|---|
| similarity | `similarity.csv`, header `model,version`, one or more models: `eos8nyg`, `eos1bgt`, `eos80al`, `eos3pu0`, `eos9b1c` or `eos8698` |
| generators | `generators.csv`, header `model,version` |
| annotators | `annotators.csv`, header `source,model,version,cutoff,direction,priority` (+ `weight,column`) |
| seeds | `smiles.csv`, header `smiles`, one SMILES per row |

`version` is a release tag of the model's GitHub repository, e.g. `v1.0.0` (its
releases are at `https://github.com/ersilia-os/<model>/releases`), or `live`: the
latest commit of its main branch, which need not be a release. A Hub model needs
one: an empty cell is an error that says where to find it. An
`rdkit` annotator has none (leave the cell empty). A model listed in two files
must have the same version in both. `fetch` checks each clone out at its
version (`live`: its main branch, pulled), and `run` refuses a clone at another
one (for `live`, one that is not on its main branch; only `fetch` can tell
whether it is the latest commit). The old one-column headers
(`similarity`, `generators`) and an `annotators.csv` without `version` are
refused. In a config file the lists are objects:
`"generators": [{"model": "eos9p57", "version": "v1.0.0"}]`, and an annotator
object takes a `version` key.

`run` needs a seed: the molecules in `smiles.csv`, or the one SMILES of
`--seed-smiles`, which overrides the file (with a note). Each seed is climbed in
turn, as if it were a run of its own, in its own folder of the run: `seed1/`,
`seed2/`, ... (zero-padded from 10 seeds on), each with its `round<n>.csv`,
`summary.csv` and, on a cluster, its `jobs/` and `jobs.json`. A stopped run
resumes in the seed it was in; the seeds before it are only read back. A SMILES
RDKit cannot read, one with a space in it (RDKit would read the rest as a
name), or the same molecule twice is refused. `init` writes no `smiles.csv`.

## Annotators

`source` is `ersilia` (`model` is a Hub id) or `rdkit` (`model` is `mw`, `qed` or
`tanimoto`). `cutoff` is a number, or `min-max` with `direction` `range`
(`higher`/`lower` take a single number, and only a guard-rail can be a range).
`priority` is `primary` (optimized, in list order for `sequential` and
`incremental`) or `guard-rail` (a constraint every winner must satisfy, no weight).

`tanimoto` is the similarity of each candidate to the seed of its climb (Morgan
fingerprints, radius 2, 2048 bits; the seed itself scores 1.0).
As a guard-rail it keeps candidates near the seed (`0.4,higher`), away from it
(`0.6,lower`) or both (`0.4-0.9,range`); as a primary, `lower` pushes the search
away from the seed. `mw`, `qed` and `tanimoto` run in the `rdkit` environment,
which `chemsampler fetch rdkit` builds (this Python's version and chemsampler's
own RDKit, so the values are the same everywhere): as a process of its own
here, as a job on a cluster. `weight` (primary, default 1) and `column` (a
multi-output Hub model) are optional.

## The project file

`chemsampler init` writes every location into the project's `chemsampler.json`,
as absolute paths:

```json
{
  "path_to_models": "/home/you/projects/my-run/input",
  "path_to_runs": "/home/you/projects/my-run/runs",
  "path_to_model_repositories": "/home/you/.chemsampler/models",
  "path_to_envs": "/home/you/.chemsampler/envs"
}
```

`path_to_models` is the folder with the input CSVs, `path_to_runs` where each run
writes its dated folder, and the last two where `fetch` puts the model clones and
their conda environments: this machine's store unless `init --models DIR
--envs DIR` named others. To change a location later, edit the file (a relative
path there is relative to the project folder); never move model or environment
folders after `fetch`.

A copied or moved project still names the old folders. Commands refuse a project
whose input or runs folder belongs to another project, so a copy never reads the
original's inputs or writes into its results: fix the two paths in the copy's
`chemsampler.json`.

## The catalog

Every Hub model is checked against
[`src/chemsampler/hub/catalog.csv`](../src/chemsampler/hub/catalog.csv): generators
must be labelled Sampling / Generation, the similarity model Sampling / Similarity
search, annotators Annotation. The catalog lists only models the Hub lists (the
Ready ones). Of the Hub's similarity models only six are accepted: `eos8nyg` (Ersilia
Reference Library), `eos1bgt` (COCONUT natural products), `eos80al` (Enamine Hit
Locator), `eos3pu0` (Enamine Liquid Stock), `eos9b1c` (Molport) and `eos8698`
(Enamine REAL sample); the others call outside services. The check is instant and works
offline.

`chemsampler catalog --update` rebuilds the catalog from the Hub (about 20 s, needs
a network connection), tells you which models are new or gone, and leaves the file
untouched if the Hub cannot be reached. Commit the file after running it.

A similarity-search model returns the seed's nearest neighbours in a reference
library: known compounds, not generated ones. Its neighbours of the current best
molecule join each round's candidates, and the `source` column shows its id.
