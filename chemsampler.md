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

## Role: guide, not decision-maker

**This skill only guides the user through the chemsampler pipeline. It never
makes a decision that belongs to the user.** Every choice that shapes a run —
which annotators, cutoffs and directions; the combination mode; the seed; any
Tanimoto constraint; which generators; the number of rounds, tolerance and
backend; whether to run at all; whether to run another round — is the user's.

What the skill does:
- **Explains** what each choice means and what it will cost or change (a
  mode's effect, a generator's measured behaviour and time, whether a cutoff
  is achievable).
- **Warns** when something looks unrealistic or is likely not to work, and
  says why.
- **Recommends**, when it has grounds to. A recommendation is labelled as
  one, comes with its reasoning and evidence (measured figures over
  documentation), and always leaves the user free to go another way.
- **Shows** every value that will be used, defaults included, before
  anything is written or run.

What the skill never does:
- Silently pick, default or infer a value the user has not seen and agreed
  to. A default is a proposal until the user accepts it.
- Drop, add or reorder annotators or generators, change a cutoff, or switch
  backend on its own, even when its own analysis says that would help. It
  says so and lets the user decide.
- Run anything — including the baseline check in step 6 and any follow-up
  round — without the user's go-ahead.
- Override a choice the user made after being warned. The single exception
  is a value that is mathematically impossible (step 2), which the skill
  declines to write; even then the user supplies the corrected value.

When unsure whether something is the user's call, treat it as theirs and ask.

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
     mathematically impossible to satisfy. **Decline to write it, explain
     why, and ask the user for a corrected value.** This is the one case
     where the skill won't follow a value as given: no run can ever satisfy
     it, so the row would only make the run silently fail or return nothing.
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

3. **Combination mode** — always ask, whatever the number of annotators (the
   CLI requires `--mode` and has no default, and with exactly one annotator
   the two modes diverge sharply). Explain both. If the user's stated goal
   clearly points to one (e.g. "get the best possible value" with a single
   annotator points to `sequential`), recommend it; the user picks:
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
   write a SMILES from memory (see `CLAUDE.md`) — and show the user what it
   resolved to so they can confirm it is the compound they meant. If no seed
   is given, only
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

6. **Baseline check & generator recommendation.** Before picking
   generators or running anything substantial, offer to get the seed's
   current score against each annotator, and do it once the user agrees (it
   is a run against the Hub, so step 10's rule applies):
   `chemsampler run --n-rounds 0` (annotators only;
   `--generators` can stay at its default, since no generator is ever
   invoked when `n_rounds=0` — the round loop simply doesn't execute, only
   the seed's own round-0 score is computed) returns immediately with that
   baseline via `hill_climb`'s existing round-0 behavior. Compare it to each
   annotator's cutoff and use the gap to inform a generator
   **recommendation**, which the user is free to override:
   - **The 6 shipped generators, in the user's two groups** (mean Tanimoto to
     input and time per seed from the shared benchmark in
     ersilia-os/ersilia#1919; times are best-case, amortised over a
     100-compound batch, and generators run sequentially, so they add up to
     ~76 s per seed):
     - *Move around chemical space*: `eos8vud` (0.15, ~7 s), `eos2401` (0.13,
       ~9 s), `eos935d` (0.63, ~4 s).
     - *Explore the near surroundings*: `eos8fma` (0.68, ~10 s), `eos4q1a`
       (0.61, ~7 s), `eos633t` (0.33, ~40 s).

     The groups are the user's choice, not a measurement: `eos935d` measures
     as near (0.63) and `eos633t` as in between (0.33, although it keeps the
     scaffold in 100% of cases). `--tanimoto-cutoff` gates on fingerprint
     Tanimoto, so when a floor is set use the measured figure, not the group
     label.
   - **Large gap, no (or a loose) Tanimoto floor** → recommend the
     generators that make bigger structural jumps (`eos8vud`, `eos2401`; the
     latter keeps only small 60-100 Da ring fragments rather than rebuilding
     from the seed's scaffold) — scaffold-preserving generators may not be
     able to close a large gap.
   - **A binding Tanimoto floor is set** (e.g. `>= 0.4`) → recommend against
     prioritizing generators likely to produce candidates below that floor,
     regardless of how large the annotator gap is — a candidate that fails
     the similarity gate is wasted no matter its score, so the floor takes
     precedence over the gap, and a slow generator that mostly gets filtered out is doubly
     wasted (`eos633t` is the slowest shipped generator and, at 0.33 mean
     Tanimoto, the likeliest to be gated out; see also the eos6ost note
     below).
   - **Low-Tanimoto generators to flag under a binding floor**:
     `eos8vud` (0.15) and `eos2401` (0.13) on the #1919 benchmark, and to a
     lesser degree `eos633t` (0.33).
   - **A recorded contradiction, `eos6ost` (no longer in the shipped
     default)**: live-tested twice (aspirin seed, two different Tanimoto
     floors), its candidates landed almost entirely at Tanimoto 0.04-0.1 to
     the seed, and it was by far the slowest generator (~15-22 minutes vs.
     well under a minute for the others combined). The #1919 benchmark
     disagrees on both counts: scaffold kept in 100% of cases, mean
     Tanimoto 0.27, ~35 s per compound. Unresolved. One unverified
     hypothesis for the Tanimoto part: aspirin's Murcko scaffold is only
     benzene, so a scaffold-preserving generator can still land far from it
     in fingerprint space. Treat `eos6ost`'s behaviour as seed-dependent
     until measured, and prefer measurement over its docstring either way.
   - **State this as a recommendation, not a fact, but weight real
     evidence over documentation when you have it.** Absent measured data,
     this is qualitative guidance from each generator's documented behavior
     (see the docstring in `src/chemsampler/models/generator.py`) — say so,
     and let the user override it. The 6 shipped generators have measured
     Tanimoto figures (#1919); any other Hub generator does not, so say so
     when recommending one. When step 11's empirical
     `source`/`tanimoto_to_seed` history is available from a prior run
     against the same or a similar seed, use it directly instead of
     guessing from docs — that's exactly how the `eos6ost` finding above
     was made.

7. **Generators.** Propose the shipped 6-generator CSV, adjusted by step 6's
   recommendation when a seed was given, and show the list. The user accepts
   it or removes/adds generators; the final list is theirs, and the skill
   never drops or adds one on its own. Whatever the user confirms is what
   step 9 writes to `generators.csv`.

8. **Remaining knobs.** Don't interrogate the user about each one, but don't
   apply any silently either: state the value that will be used, so it is
   visible and easy to change (step 10 shows them again in the final
   command):
   - `--n-rounds` (default 5)
   - `--tolerance` (default 0.0)
   - `--backend` (default `run_sh`, which needs the model already fetched
     locally and conda-packed; `ersilia` is the fallback). If `run_sh` fails
     because a model isn't available locally, tell the user and let them
     choose whether to fetch it or fall back to `ersilia`; never switch on
     its own.

9. **Output directory.** Propose a fresh, dedicated `--output-dir` per
   invocation (e.g. timestamped) — `chemsampler run` never clears stale
   files in a reused directory, so reusing one across runs can leave old
   `round4.csv`/`round5.csv` behind from an earlier, longer run. If the user
   wants to reuse a directory, explain that risk and let them decide. Write the
   `annotators.csv`/`generators.csv` built above into that same directory
   so each run is self-contained and reproducible.

10. **Confirm before running.** Show the exact `chemsampler run` command,
    with every value that will be used (defaults included), before executing
    it — a run can be slow and hits the real Ersilia Model Hub. Run it only
    after the user's explicit go-ahead, and ask again if anything changed.
    The same applies to the baseline check in step 6 and to any follow-up
    round from step 11.

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
      round" vs. "cutoff reached / no improvement — stop here"). This is a
      suggestion; whether to run another round is the user's call.

    Point to the per-round CSVs for full detail rather than reproducing
    them wholesale in chat. Also close the loop on step 6's generator
    guidance empirically: each round's CSV already has a `source` column
    (which generator(s) produced each candidate) and, if a seed was given,
    `tanimoto_to_seed`. If one generator's candidates are consistently
    filtered out by the Tanimoto gate, or it contributed nothing at all,
    tell the user and suggest dropping it from `generators.csv` for a
    re-run rather than silently carrying dead weight; the user decides.

## Open questions (not yet resolved — do not guess)

- **Patience-based stopping.** A possible alternative/complement to
  `cutoff`: "stop if the value hasn't improved by X in the last N rounds."
  Not designed; see session notes.
- **Final skill location.** `.claude/skills/` (machine-local, gitignored)
  vs. the shared `ersilia-skills` repo — not decided.
