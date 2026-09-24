# chemsampler skill (draft placeholder)

> **Status: draft.** This is a first-draft description of a Claude Code skill
> that coordinates `chemsampler` conversationally. It is not yet wired up as
> an active skill — it needs to be copied into an actual skill location
> (`.claude/skills/` for a machine-local skill, or the shared
> `ersilia-skills` repo) before it can be invoked. Written here so it lives
> in version control and is easy to review/iterate on.

## Purpose

Help a chemist run `chemsampler run` without needing to hand-write the
annotators/generators CSVs or remember the CLI flags. The skill's job is to
ask a small number of questions, build the input CSVs, invoke the CLI, and
report the result — not to reimplement any optimization logic itself.

**The skill drives chemsampler through the `chemsampler run` CLI, not the
Python API.** This matches how every other Ersilia skill shells out to its
tool's own CLI rather than importing its Python classes directly.

## Conversation flow

1. **Annotators first.** Start with "what are you optimizing for?", not "do
   you have a starting molecule?" — a chemist typically knows their target
   property before they know a seed. For each annotator the user names,
   collect:
   - `annotator_id` — an Ersilia Hub model id, or `qed`.
   - `cutoff` — the threshold value.
   - `direction` — `higher` or `lower`.

   All annotators are treated uniformly; there is no directing/controlling
   distinction. Do not write `annotators.csv` yet — run the feasibility
   check below first.

2. **Feasibility & relevance check.** For every cutoff collected above (and
   for `--tanimoto-cutoff` in step 5, once given), judge it two distinct
   ways before it ever reaches `annotators.csv` or the CLI:

   - **Feasible?** Can the model's output possibly take this value at all?
     Some outputs have a hard, known range regardless of which specific
     model produces them — anything framed as a probability (an
     inhibition/binding probability, a MAIP-style score, etc.) is bounded
     to `[0, 1]`; Tanimoto similarity is likewise always `[0, 1]`. A cutoff
     outside a value's hard range (e.g. "inhibition probability below
     -0.1", or a Tanimoto cutoff of 1.2) is not just unlikely, it is
     mathematically impossible to satisfy. **Refuse and explain why** —
     don't write a CSV row that can never be satisfied and let the run
     silently fail or return nothing.
   - **Relevant?** Even when a value is technically achievable, is it
     realistic for what the Ersilia Model Hub's models were built and
     validated on? Open-ended physicochemical descriptors (molecular
     weight, logP, etc.) have no hard bound, so "MW < 50 Da" or "MW >
     10,000 Da" is technically satisfiable, but the Hub's generators were
     validated on drug-like chemical space (roughly 200-450 Da for MW, per
     the project's own generator-validation notes) — a cutoff far outside
     that band is unlikely to be meaningful. **Warn explicitly** (name the
     applicability-domain concern) and get the user to confirm before
     proceeding, rather than silently accepting it.
   - **How to judge this**: there's no code-level registry of per-model
     output ranges in chemsampler today, so this is conversational
     judgment, not a library check. Recognize known output *types*
     (probability-like outputs are always `[0, 1]` no matter which eos
     code produces them); apply general domain knowledge about
     physicochemical descriptor ranges for drug-like chemistry; and if
     genuinely unsure what an unfamiliar Hub model's output means or its
     plausible range, say so and ask the user rather than silently trusting
     an extreme cutoff.

3. **Combination mode**, only if more than one annotator was given (with
   exactly one annotator the two modes diverge sharply, so still ask
   explicitly rather than defaulting):
   - `sequential` — annotators are optimized one at a time, in priority
     order (the order the user listed them). Each finished stage's achieved
     value becomes a hard floor for every later stage. Reproduces plain
     maximization/minimization when there's a single annotator.
   - `joint` — all annotators optimized simultaneously; a candidate's score
     is how many cutoffs it satisfies. With a single annotator this is a
     threshold + arbitrary pick among passers, not maximization — mention
     this explicitly if the user has only one annotator and seems to expect
     "get the best possible value."

4. **Seed compound (optional).** Ask for a SMILES. If the user names a
   compound instead of pasting a SMILES, look it up from an authoritative
   source (PubChem) and confirm by **InChIKey** before using it — never
   write a SMILES from memory (see `CLAUDE.md`). If no seed is given, only
   generators that don't require one will produce candidates in round 1,
   and steps 5-6 below don't apply (no seed to be similar to, or to score a
   baseline against).

5. **Tanimoto constraint (optional, only if a seed was given).** Ask
   whether to bound similarity to the seed at all, and if so
   `--tanimoto-cutoff` + `--tanimoto-direction` (`higher` = stay at least
   this similar; `lower` = novelty-seeking, stay at most this similar). Run
   it through the feasibility check from step 2 — it's always `[0, 1]`.
   Collected here, before generator selection, because it directly
   constrains which generators make sense to use next.

6. **Baseline check & generator selection guidance.** Before picking
   generators or running anything substantial, get the seed's current score
   against each annotator: `chemsampler run --n-rounds 0` (annotators only;
   `--generators` can stay at its default, since no generator is ever
   invoked when `n_rounds=0` — the round loop simply doesn't execute, only
   the seed's own round-0 score is computed) returns immediately with that
   baseline via `hill_climb`'s existing round-0 behavior. Compare it to each
   annotator's cutoff and use the gap to inform generator choice:
   - **Large gap, no (or a loose) Tanimoto floor** → lean toward generators
     documented to make bigger structural jumps (e.g. `eos2401`, which keeps
     only small 60-100 Da ring fragments rather than rebuilding from the
     seed's scaffold) — the scaffold-preserving generators (`eos9taz`,
     `eos6ost`) may not be able to close a large gap.
   - **A binding Tanimoto floor is set** (e.g. `>= 0.5`) → don't prioritize
     generators likely to produce candidates below that floor, regardless of
     how large the annotator gap is — a candidate that fails the similarity
     gate is wasted no matter its score, so the floor takes precedence over
     the gap. Prefer the scaffold-preserving generators here.
   - **State this as a recommendation, not a fact.** This is qualitative
     guidance from each generator's documented behavior (see the docstring
     in `src/chemsampler/models/generator.py`), not measured Tanimoto
     statistics — chemsampler doesn't track empirical similarity
     distributions per generator. Say so, and let the user override it.

7. **Generators (optional).** Default to the shipped 3-generator CSV,
   informed by step 6's guidance when a seed was given; ask only if the
   user wants to customize further. If customizing, write a fresh
   `generators.csv`.

8. **Remaining knobs**, offered with sensible defaults rather than asked
   outright unless the user wants to tune them:
   - `--n-rounds` (default 5)
   - `--tolerance` (default 0.0)
   - `--backend` (default `ersilia`; `run_sh` is an alternative execution
     path, not something to offer unless asked)

9. **Output directory.** Use a fresh, dedicated `--output-dir` per
   invocation (e.g. timestamped) — `chemsampler run` never clears stale
   files in a reused directory, so reusing one across runs can leave old
   `round4.csv`/`round5.csv` behind from an earlier, longer run. Write the
   `annotators.csv`/`generators.csv` built above into that same directory
   so each run is self-contained and reproducible.

10. **Confirm before running.** Show the exact `chemsampler run` command
    before executing it — a run can be slow and hits the real Ersilia Model
    Hub. Don't invoke it silently.

11. **Report the result and check generator performance.** Don't just dump
    `summary.csv`/`round<n>.csv` — interpret them into a short, structured
    summary, e.g.:
    - Total number of candidates generated (across all generators, that round).
    - Per generator: how many candidates it contributed; for any generator
      that contributed 0 (chemsampler already logs a warning for this),
      say so explicitly and explain why if known (e.g. a generator that
      only accepts certain scaffold types, or the seed being outside its
      working range) rather than silently omitting it.
    - The best score achieved (and which annotator it's against, in
      `mode="sequential"`), plus that candidate's `tanimoto_to_seed` if a
      seed was given.
    - Which generator took the most wall-clock time, if that's visible
      (e.g. from timing chemsampler itself logs, or from run duration) —
      useful for the user to know which model is the bottleneck when
      iterating.
    - A concrete recommendation: given the best candidate found, does it
      make sense to run a second round using it as the new seed? (e.g.
      "score improved but is still far from the cutoff — worth another
      round" vs. "cutoff reached / no improvement — stop here").

    Point to the per-round CSVs for full detail rather than reproducing
    them wholesale in chat. Also close the loop on step 6's generator
    guidance empirically: each round's CSV already has a `source` column
    (which generator(s) produced each candidate) and, if a seed was given,
    `tanimoto_to_seed`. If one generator's candidates are consistently
    filtered out by the Tanimoto gate, or it contributed nothing at all,
    tell the user and suggest dropping it from `generators.csv` for a
    re-run rather than silently carrying dead weight.

## Open questions (not yet resolved — do not guess)

- **Patience-based stopping.** A possible alternative/complement to
  `cutoff`: "stop if the value hasn't improved by X in the last N rounds."
  Not designed; see session notes.
- **Final skill location.** `.claude/skills/` (machine-local, gitignored)
  vs. the shared `ersilia-skills` repo — not decided.
