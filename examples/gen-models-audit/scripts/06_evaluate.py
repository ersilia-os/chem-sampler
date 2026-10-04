"""
Evaluate one model's outputs over the ChEMBL seed splits.

Reads --path-to-output/<model-id>/csv/split_<NNN>.csv (written by
04_run_split.py) against results/ChEMBL_splits/split_<NNN>.csv, and reports
over every compound found (a row is one compound, a cell one output slot):

- % compounds with >= 1 null, with >= --null-threshold nulls (default 10),
  and with every slot null (100%). A null is an empty cell (or "nan" /
  "none" / "null").
- Max, P95 and mean Tanimoto to the compound's own input: for each compound
  the maximum, the 95th percentile (linear interpolation) and the mean over
  its valid outputs; then the mean +- sample standard deviation (ddof=1) of
  each across the compounds that have at least one valid output. Morgan
  fingerprints, radius 2, 2048 bits.
- Retention of the input, per compound as a share of its valid outputs, then
  mean +- std over the compounds that have outputs: % with the same Murcko
  scaffold as the input, % with the same generic scaffold (every atom a
  carbon, every bond single), % containing the whole input as a substructure
  (all-structure), and % with more heavy atoms than the input. Compounds whose
  input has no ring are left out of the two scaffold figures. The same three
  levels are also read together as a breakdown that adds up to 100%: each
  output counts once, at the most specific level it keeps (the whole input,
  else the Murcko scaffold, else the generic scaffold, else none).
- Quality checks over all non-null outputs: invalid SMILES, duplicates
  (canonical isomeric SMILES, within a compound), input echoes (same
  structure as the input ignoring stereochemistry), multi-component ("."),
  atom-map numbers and dummy atoms ("*").
- Wall time from the run provenance in logs/*.json.

If split_<NNN>_repeated.csv files exist, the same split is also evaluated
from each run, next to how many molecules the two runs share (identical
canonical SMILES) per compound and whether the same compounds come back
empty. Writes results/<model-id>/analysis.json and md/analysis.md and prints
the markdown.

Needs rdkit and numpy. Fully self-contained: no imports from chemsampler.
"""

import argparse
import csv
import glob
import json
import os
import re

import numpy as np
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
AUDIT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
DEFAULT_SPLITS_DIR = os.path.join(AUDIT_DIR, "results", "ChEMBL_splits")
DEFAULT_OUTPUT_DIR = os.path.join(AUDIT_DIR, "results")

NULL_STRINGS = {"", "nan", "none", "null"}
ATOM_MAP_RE = re.compile(r":\d+\]")
FP_RADIUS = 2
FP_BITS = 2048
FP_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=FP_RADIUS, fpSize=FP_BITS)


def read_rows(path: str) -> list[list[str]]:
    """Data rows of a CSV file (header dropped)."""
    with open(path, newline="") as f:
        return list(csv.reader(f))[1:]


def is_null(cell: str) -> bool:
    return cell.strip().lower() in NULL_STRINGS


def canonical(mol, isomeric: bool = True) -> str:
    return Chem.MolToSmiles(mol, isomericSmiles=isomeric)


def generic_scaffold_smiles(scaffold) -> str | None:
    """SMILES of a Murcko scaffold made generic (all carbons, single bonds), or
    None when RDKit cannot build it (e.g. a hypervalent P or S becomes a
    five- or six-valent carbon)."""
    try:
        return Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(scaffold))
    except Exception:  # noqa: BLE001
        return None


def analyse_compound(seed: str, cells: list[str]) -> dict:
    """Per-compound counts and Tanimoto statistics for one row of outputs."""
    seed_mol = Chem.MolFromSmiles(seed)
    seed_flat = canonical(seed_mol, isomeric=False)
    seed_fp = FP_GEN.GetFingerprint(seed_mol)
    seed_scaffold = MurckoScaffold.GetScaffoldForMol(seed_mol)
    has_scaffold = seed_scaffold.GetNumAtoms() > 0
    seed_scaffold_smiles = Chem.MolToSmiles(seed_scaffold) if has_scaffold else None
    seed_generic_smiles = (
        generic_scaffold_smiles(seed_scaffold) if has_scaffold else None
    )
    has_generic = seed_generic_smiles is not None
    seed_heavy = seed_mol.GetNumHeavyAtoms()
    kept = {"scaffold": 0, "generic_scaffold": 0, "all_structure": 0, "heavier": 0}
    profile = {"whole": 0, "murcko": 0, "generic": 0, "none": 0}

    outputs = [c.strip() for c in cells if not is_null(c)]
    seen: set[str] = set()
    result = {
        "n_slots": len(cells),
        "n_null": len(cells) - len(outputs),
        "invalid": 0,
        "duplicates": 0,
        "echoes": 0,
        "multi_component": 0,
        "atom_map": 0,
        "dummy_atom": 0,
        "canonical": set(),
    }
    fps = []
    for smi in outputs:
        result["multi_component"] += "." in smi
        result["atom_map"] += bool(ATOM_MAP_RE.search(smi))
        result["dummy_atom"] += "*" in smi
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            result["invalid"] += 1
            continue
        key = canonical(mol)
        if key in seen:
            result["duplicates"] += 1
        seen.add(key)
        result["echoes"] += canonical(mol, isomeric=False) == seed_flat
        fps.append(FP_GEN.GetFingerprint(mol))
        same_scaffold = same_generic = False
        if has_scaffold:
            scaffold = MurckoScaffold.GetScaffoldForMol(mol)
            same_scaffold = Chem.MolToSmiles(scaffold) == seed_scaffold_smiles
            if has_generic:
                same_generic = generic_scaffold_smiles(scaffold) == seed_generic_smiles
        whole = mol.HasSubstructMatch(seed_mol)
        kept["scaffold"] += same_scaffold
        kept["generic_scaffold"] += same_generic
        kept["all_structure"] += whole
        kept["heavier"] += mol.GetNumHeavyAtoms() > seed_heavy
        # the most specific thing this output keeps from its input
        if whole:
            profile["whole"] += 1
        elif same_scaffold:
            profile["murcko"] += 1
        elif same_generic:
            profile["generic"] += 1
        else:
            profile["none"] += 1
    result["canonical"] = seen
    available = {
        "scaffold": has_scaffold,
        "generic_scaffold": has_generic,
        "all_structure": True,
        "heavier": True,
    }
    result["retention"] = (
        {
            key: count / len(fps) if available[key] else None
            for key, count in kept.items()
        }
        if fps
        else None
    )
    result["profile"] = (
        {level: count / len(fps) for level, count in profile.items()} if fps else None
    )

    if fps:
        sims = np.array(DataStructs.BulkTanimotoSimilarity(seed_fp, fps))
        result["tanimoto"] = {
            "max": float(sims.max()),
            "p95": float(np.percentile(sims, 95)),
            "mean": float(sims.mean()),
        }
    else:
        result["tanimoto"] = None
    return result


def mean_std(values: list[float]) -> dict:
    """Mean and sample standard deviation (ddof=1; 0 for a single value)."""
    arr = np.array(values, dtype=float)
    if arr.size == 0:
        return {"mean": None, "std": None, "n": 0}
    std = float(arr.std(ddof=1)) if arr.size > 1 else 0.0
    return {"mean": float(arr.mean()), "std": std, "n": int(arr.size)}


def summarise(compounds: list[dict], null_threshold: int) -> dict:
    """Aggregate per-compound results into the reported metrics."""
    n = len(compounds)
    if n == 0:
        return {"n_compounds": 0}

    def pct(count: int) -> float:
        return 100.0 * count / n

    def total(key: str) -> int:
        return sum(c[key] for c in compounds)

    n_outputs = sum(c["n_slots"] - c["n_null"] for c in compounds)
    tan = [c["tanimoto"] for c in compounds if c["tanimoto"]]

    def retention(key: str) -> dict:
        shares = [
            100.0 * c["retention"][key]
            for c in compounds
            if c["retention"] and c["retention"][key] is not None
        ]
        return mean_std(shares)

    return {
        "n_compounds": n,
        "n_outputs": n_outputs,
        "pct_with_any_null": pct(sum(c["n_null"] >= 1 for c in compounds)),
        f"pct_with_{null_threshold}_or_more_nulls": pct(
            sum(c["n_null"] >= null_threshold for c in compounds)
        ),
        "pct_all_null": pct(sum(c["n_null"] == c["n_slots"] for c in compounds)),
        "max_tanimoto": mean_std([t["max"] for t in tan]),
        "p95_tanimoto": mean_std([t["p95"] for t in tan]),
        "mean_tanimoto": mean_std([t["mean"] for t in tan]),
        "pct_scaffold": retention("scaffold"),
        "pct_generic_scaffold": retention("generic_scaffold"),
        "pct_all_structure": retention("all_structure"),
        "pct_heavier": retention("heavier"),
        "retention_profile": {
            level: mean_std(
                [100.0 * c["profile"][level] for c in compounds if c["profile"]]
            )
            for level in ("whole", "murcko", "generic", "none")
        },
        "mean_outputs_per_compound": n_outputs / n,
        "invalid": total("invalid"),
        "duplicates": total("duplicates"),
        "echoes": total("echoes"),
        "multi_component": total("multi_component"),
        "atom_map": total("atom_map"),
        "dummy_atom": total("dummy_atom"),
    }


def load_split_seeds(splits_dir: str, name: str) -> list[str]:
    return [row[0] for row in read_rows(os.path.join(splits_dir, f"{name}.csv"))]


def evaluate_files(paths: list[str], splits_dir: str) -> list[dict]:
    """Per-compound results for every output file in `paths`."""
    compounds = []
    for path in paths:
        name = os.path.basename(path)[: -len(".csv")].replace("_repeated", "")
        seeds = load_split_seeds(splits_dir, name)
        rows = read_rows(path)
        if len(rows) != len(seeds):
            raise SystemExit(f"{path}: {len(rows)} rows for {len(seeds)} inputs")
        compounds += [analyse_compound(s, r) for s, r in zip(seeds, rows)]
    return compounds


def wall_times(logs_dir: str) -> dict:
    """Wall-time summary from the per-run provenance files."""
    runs = []
    for path in sorted(glob.glob(os.path.join(logs_dir, "split_*.json"))):
        with open(path) as f:
            runs.append(json.load(f))
    ok = [r for r in runs if r["exit_code"] == 0 and not r["repeated"]]
    if not ok:
        return {}
    secs = [r["wall_seconds"] for r in ok]
    return {
        "n_runs": len(ok),
        "failed_runs": sum(r["exit_code"] != 0 for r in runs),
        "mean_seconds_per_split": float(np.mean(secs)),
        "seconds_per_100_compounds": float(np.mean(secs)) * 100 / ok[0]["n_inputs"],
        "hosts": sorted({r["host"].split(".")[0] for r in ok}),
        "cpus": sorted({r["cpus"] for r in ok}),
    }


def compare_runs(first: list[dict], second: list[dict]) -> dict:
    """How much two runs of the same split agree, compound by compound."""
    shared, jaccard, same_empty = [], [], 0
    for a, b in zip(first, second):
        common = a["canonical"] & b["canonical"]
        union = a["canonical"] | b["canonical"]
        shared.append(len(common))
        if union:
            jaccard.append(len(common) / len(union))
        same_empty += (a["n_null"] == a["n_slots"]) == (b["n_null"] == b["n_slots"])
    return {
        "n_compounds": len(first),
        "shared_molecules_per_compound": mean_std([float(x) for x in shared]),
        "jaccard": mean_std(jaccard),
        "same_empty_status": same_empty,
    }


def fmt(stat: dict, digits: int = 3) -> str:
    if stat["mean"] is None:
        return "n/a"
    return f"{stat['mean']:.{digits}f} ± {stat['std']:.{digits}f} (n={stat['n']})"


def metric_rows(summary: dict, threshold: int) -> list[tuple[str, str]]:
    return [
        ("Compounds", str(summary["n_compounds"])),
        ("% compounds with >= 1 null", f"{summary['pct_with_any_null']:.1f}"),
        (
            f"% compounds with >= {threshold} nulls",
            f"{summary[f'pct_with_{threshold}_or_more_nulls']:.1f}",
        ),
        ("% compounds with all slots null (100%)", f"{summary['pct_all_null']:.1f}"),
        ("Max Tanimoto (mean ± std)", fmt(summary["max_tanimoto"])),
        ("P95 Tanimoto (mean ± std)", fmt(summary["p95_tanimoto"])),
        ("Average Tanimoto (mean ± std)", fmt(summary["mean_tanimoto"])),
        (
            "% outputs with the input's scaffold (mean ± std)",
            fmt(summary["pct_scaffold"], 1),
        ),
        (
            "% with the input's generic scaffold",
            fmt(summary["pct_generic_scaffold"], 1),
        ),
        (
            "% containing the whole input (all-structure)",
            fmt(summary["pct_all_structure"], 1),
        ),
        ("% heavier than the input", fmt(summary["pct_heavier"], 1)),
        (
            "Retention profile, % of outputs: whole input / Murcko scaffold / generic scaffold / none",
            " / ".join(
                f"{summary['retention_profile'][k]['mean']:.1f}"
                for k in ("whole", "murcko", "generic", "none")
            ),
        ),
        ("Outputs per compound (mean)", f"{summary['mean_outputs_per_compound']:.1f}"),
        (
            "Invalid / duplicate / input echo",
            f"{summary['invalid']} / {summary['duplicates']} / {summary['echoes']}",
        ),
        (
            "Multi-component / atom-map / dummy atom",
            f"{summary['multi_component']} / {summary['atom_map']} / {summary['dummy_atom']}",
        ),
    ]


def markdown(
    model: str, summary: dict, times: dict, repeat: dict | None, threshold: int
) -> str:
    lines = [f"# {model}: evaluation", "", "| Metric | Value |", "|---|---|"]
    lines += [f"| {k} | {v} |" for k, v in metric_rows(summary, threshold)]
    if times:
        lines.append(
            f"| Time per 100 compounds | {times['seconds_per_100_compounds'] / 60:.1f} min "
            f"({times['n_runs']} runs on {', '.join(times['hosts'])}, "
            f"{'/'.join(map(str, times['cpus']))} CPUs) |"
        )
    if repeat:
        lines += [
            "",
            f"## Repeated split ({repeat['split']})",
            "",
            "| Metric | Run 1 | Run 2 |",
            "|---|---|---|",
        ]
        rows_a = dict(metric_rows(repeat["run1"], threshold))
        rows_b = dict(metric_rows(repeat["run2"], threshold))
        lines += [f"| {k} | {rows_a[k]} | {rows_b[k]} |" for k in rows_a]
        cmp = repeat["comparison"]
        shared = fmt(cmp["shared_molecules_per_compound"], 1)
        text = (
            f"Molecules shared by both runs, per compound: {shared}; "
            f"Jaccard: {fmt(cmp['jaccard'])}; compounds with the same "
            f"empty/non-empty status in both runs: "
            f"{cmp['same_empty_status']}/{cmp['n_compounds']}."
        )
        lines += ["", text]
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", required=True, help="Model identifier, e.g. eos2401")
    parser.add_argument("--null-threshold", type=int, default=10)
    parser.add_argument("--path-to-splits", default=DEFAULT_SPLITS_DIR)
    parser.add_argument("--path-to-output", default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()

    model_out = os.path.join(args.path_to_output, args.model)
    csv_dir = os.path.join(model_out, "csv")
    outputs = sorted(
        p
        for p in glob.glob(os.path.join(csv_dir, "split_*.csv"))
        if not p.endswith("_repeated.csv")
    )
    if not outputs:
        raise SystemExit(f"no split_*.csv in {csv_dir}")

    compounds = evaluate_files(outputs, args.path_to_splits)
    summary = summarise(compounds, args.null_threshold)
    times = wall_times(os.path.join(model_out, "logs"))

    repeat = None
    repeated = sorted(glob.glob(os.path.join(csv_dir, "split_*_repeated.csv")))
    if repeated:
        path2 = repeated[0]
        name = os.path.basename(path2)[: -len("_repeated.csv")]
        path1 = os.path.join(csv_dir, f"{name}.csv")
        if os.path.exists(path1):
            run1 = evaluate_files([path1], args.path_to_splits)
            run2 = evaluate_files([path2], args.path_to_splits)
            repeat = {
                "split": name,
                "run1": summarise(run1, args.null_threshold),
                "run2": summarise(run2, args.null_threshold),
                "comparison": compare_runs(run1, run2),
            }

    report = {
        "model": args.model,
        "splits_evaluated": len(outputs),
        "summary": summary,
        "wall_time": times,
        "repeated_split": repeat,
    }
    with open(os.path.join(model_out, "analysis.json"), "w") as f:
        json.dump(report, f, indent=2)
        f.write("\n")
    text = markdown(args.model, summary, times, repeat, args.null_threshold)
    os.makedirs(os.path.join(model_out, "md"), exist_ok=True)
    with open(os.path.join(model_out, "md", "analysis.md"), "w") as f:
        f.write(text)
    print(f"{len(outputs)} split(s) evaluated")
    print(text)


if __name__ == "__main__":
    main()
