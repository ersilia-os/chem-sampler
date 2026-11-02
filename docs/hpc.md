# Running on a Slurm cluster

`chemsampler run --hpc NAME` runs every model of a run as a Slurm job. Without
`--hpc`, models run on the machine you are on, one after another.

## What a round does

1. The similarity model and the generative models start together, one job each,
   on the current best molecule.
2. When all of them are done, their candidates are put together (as a local run
   does).
3. Every annotator starts, one job each: Hub models and the RDKit ones (`mw`,
   `qed`).
4. When all of them are done, the round is ranked and its best molecule seeds the
   next round.

The seed itself is scored first, as round 0. The terminal shows one line whenever
the jobs of a step change state:

```
Round 1 · sampling: 2 job(s): 1 PD, 1 R, 0 done (3 min)
Round 1 · sampling: 2 job(s) done (5 min)
```

Each step waits in the cluster's queue, so a run is worth sending when its models
are slow; fast models can finish sooner locally.

## The profile

`hpc.json` maps profile names to clusters. It lives in the store, chemsampler's one
folder (`$CHEMSAMPLER_HOME`, or `~/.chemsampler`), beside `models/` and `envs/`:
it is yours, not part of a project or of chemsampler. Every key is optional.

| Key | What it does |
|---|---|
| `ssh` | Host (or `~/.ssh/config` alias) of a login node with Slurm. Leave it out when chemsampler runs where `sbatch` does. Keys only: it runs with `BatchMode=yes`. |
| `partition` | `--partition` for every job. |
| `nodes` | List of nodes, passed as `--nodelist` with `--nodes=1`, so that a job may run on any of them. |
| `cpus` | `--cpus-per-task` (default 1). The thread limits of the job (`OMP_NUM_THREADS` and others) follow it. |
| `mem`, `time` | `--mem` and `--time`, as sbatch takes them (`"8G"`, `"01:00:00"`). |
| `poll_seconds` | How often the jobs are checked (default 60). |
| `models` | Resources per model, e.g. `{"eos8nyg": {"mem": "16G"}, "qed": {"cpus": 1}}`. |

Several profiles can live in the same file; pick one per run with `--hpc NAME`.

## Before the first job

`run --hpc` checks and lists every problem before submitting anything:

- every model has the environment `chemsampler fetch` built, and the `rdkit`
  environment is there when RDKit annotators are used;
- Slurm's commands answer, and the partition and nodes exist;
- the run's folders, the models and their environments are visible from the
  cluster.

The cluster has to see the same paths: put the project and the store on a shared
filesystem (`export CHEMSAMPLER_HOME=/shared/.chemsampler`, or one project's
folders with `init --models/--envs`). The nodes need neither conda nor chemsampler. A
model runs its own `run.sh` with its environment first on `PATH`; an RDKit
annotator runs a copy of chemsampler's scoring file with the `rdkit` environment,
which has chemsampler's own RDKit version, so scores equal a local run's.

## What a run leaves in its folder

Besides the usual `summary.csv`, `round<n>.csv` and `inputs.json` (whose `run`
entry names the profile), a cluster run writes:

- `jobs.json`: every job, with its step, model, Slurm id, state, node, submit,
  start and end times, run time and exit code;
- `jobs/<NN>-round<r>-<step>/<model>/`: the job's `input.csv`, `job.sh`,
  `output.csv`, and Slurm's `slurm.out` and `slurm.err`.

Jobs are named `<model>-r<round>` (e.g. `eos9p57-r1`), so `squeue -u $USER` shows
which model and round each one is.

## When things go wrong

Every round is saved as soon as it is ranked, so a stopped run is never lost:
`chemsampler resume` continues it (see below).

- **A job fails** (an error, a time limit, out of memory): the run stops. The other
  jobs of the step are cancelled, the error names the model and the end of its
  `slurm.err`, and the finished rounds are written.
- **Ctrl-C** (or SIGTERM): every job of the run still on the cluster is cancelled
  (`scancel`), and the finished rounds are written.
- **Preemption**: jobs are submitted with `--requeue`. A preempted job goes back
  to the queue and starts over; its output is only renamed into place when it is
  complete, so a half-written file is never read.
- **The terminal closes**: chemsampler is the one waiting for the jobs, so the run
  stops with it and its jobs carry on alone. Start long runs inside `tmux` or
  `screen`. If it happens anyway, `resume` cancels what the run left on the
  cluster before going on.

## Resuming a run

`chemsampler resume [RUN]` (the latest run by default) continues a run with the
inputs and options it started with, including its profile:

- The rounds it finished are read back from their `round<n>.csv`; no job is sent
  for them, and every choice that followed them (winners, stages, steps) comes out
  as before.
- The round it was in is run again whole. Its old job folders move to
  `jobs/stopped/attempt<n>/`, and `jobs.json` keeps every attempt's jobs.
- Jobs an earlier attempt left on the cluster are cancelled first.
- `run.lock` in the run folder says who runs it. `resume` refuses a run that is
  still going (a live process, or a Slurm job still queued or running) and says
  how to stop it.
