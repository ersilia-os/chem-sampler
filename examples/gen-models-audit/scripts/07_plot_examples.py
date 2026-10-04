"""
Draw a model's example compounds next to their most similar generated outputs.

For each compound in the model's examples/run_input.csv (3 for most models)
one row of a grid image is drawn: column 1 is the input molecule, labelled
"INPUT <row number>", and columns 2-6 are the 5 distinct valid outputs with the highest
Tanimoto similarity (TS) to it (Morgan radius 2, 2048 bits, as in
06_evaluate.py), most similar first, each labelled "TS=<value> (top-<rank>)".
A compound with fewer than 5 valid outputs leaves the remaining cells blank,
labelled "no output" in grey.

Highlighting follows the model's Type in gen-models-master.csv (override with
--highlight):
- Scaffold-based: the input's Murcko scaffold, in red, in the input and in
  every output that still contains it;
- Fragment-seeded growth / Fragment recombination: the fragment of the
  input that is kept in the outputs, in blue. For each output it is the
  complete rings of its largest common substructure with the input (the whole
  common substructure if that has no ring, at least 3 atoms), so chain atoms
  that merely happen to match are left out; the input shows the fragments
  kept by any of its outputs, nothing else;
- any other type: no highlighting.

Outputs come from the model's shipped examples/run_output.csv by default, so
the picture shows what the Hub ships and needs no compute; pass --outputs to
draw a fresh run instead (a CSV with one row per example compound).
The image is written to --path-to-output/<model-id>/png/examples_grid.png.

A second figure, --path-to-output/<model-id>/png/chemical_space.png, maps the same
compounds in chemical space: the t-SNE coordinates of the Ersilia model eos1klk
(a projector trained on the Ersilia reference library) are computed for 20,000
random compounds of the reference library (lazy-chemvis data/smiles_100k.csv,
seed 42), drawn in gray as the background, and for the example's seed compounds
(black, numbered) and up to 100 outputs per seed (crimson). There is one square
panel per seed, side by side, each with the background, its seed and its own
outputs; all layers are plain points, under a centred title (model id | slug |
Type), with no legend or grid, and dot area grows with the local density of
the layer within the panel. It is drawn with stylia (print format, article style
and colours). The
background is computed once and cached in --path-to-cache/eos1klk/. Only the t-SNE part of
eos1klk is run: for a new molecule its t-SNE coordinates are a Morgan
fingerprint (radius 2, 2048 bits) fed to one XGBoost model, so the descriptors,
PCA, TMAP and UMAP steps of the full model are skipped; the coordinates are
identical to the ones the full model returns and about 300x faster. eos1klk
must be fetched and its env built first: 02_fetch_generative_models.py
--models eos1klk, then 03_build_cpu_envs.py --models eos1klk. --plots chooses
grid, space or both.

Needs rdkit, Pillow and numpy, plus stylia (with matplotlib) for the map.
Fully self-contained: no imports
from chemsampler.
"""

import argparse
import csv
import os
import random
import subprocess
import tempfile
import urllib.request
from unittest import mock

import numpy as np
from PIL import ImageDraw, ImageFont
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Draw, rdFingerprintGenerator, rdFMCS
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
AUDIT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
DEFAULT_MODELS_DIR = os.path.join(AUDIT_DIR, "models")
DEFAULT_OUTPUT_DIR = os.path.join(AUDIT_DIR, "results")
DEFAULT_TABLE = os.path.join(AUDIT_DIR, "gen-models-master.csv")
DEFAULT_ENVS_DIR = os.path.join(AUDIT_DIR, "envs_cpu")
DEFAULT_CACHE_DIR = os.path.join(AUDIT_DIR, "cache")

PROJECTOR = "eos1klk"
REFERENCE_URL = "https://raw.githubusercontent.com/ersilia-os/lazy-chemvis/main/data/smiles_100k.csv"
N_REFERENCE = 20000
REFERENCE_SEED = 42
MAX_OUTPUTS_PER_SEED = 100
DENSITY_BINS = 100
DENSITY_SIGMA = 2.0  # in bins
# Map dot areas as multiples of stylia.MARKERSIZE; the two values are the sparse
# and the dense end of the density range
REFERENCE_SIZES = (0.1, 0.6)
OUTPUT_SIZES = (0.4, 2.0)
SEED_SIZE = 3.0
MAP_HEIGHT = 0.42  # figure height as a fraction of the width: one row of squares
TITLE_GAP_INCHES = 0.08  # between the map title and the panels
TITLE_FONT_FACTOR = 1.5  # map title size as a multiple of stylia.FONTSIZE_BIG

N_SIMILAR = 5
CELL_PIXELS = 220
FP_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
NULL_STRINGS = {"", "nan", "none", "null"}

TYPE_TO_MODE = {
    "Scaffold-based": "scaffold",
    "Fragment-seeded growth": "fragment",
    "Fragment recombination": "fragment",
}
MODE_COLOUR = {"scaffold": (0.9, 0.0, 0.0), "fragment": (0.0, 0.35, 1.0)}
MIN_FRAGMENT_ATOMS = 3
HIGHLIGHT_BOND_WIDTH = 16  # RDKit default 8
HIGHLIGHT_RADIUS = 0.4  # RDKit default 0.3
LEGEND_FONT_SIZE = 22  # RDKit default 16
MCS_TIMEOUT_SECONDS = 2


def read_rows(path: str) -> list[list[str]]:
    """Data rows of a CSV file (header dropped)."""
    with open(path, newline="") as f:
        return list(csv.reader(f))[1:]


def mode_from_table(model: str, table_path: str) -> str:
    """Highlight mode for `model` from its Type in the master table."""
    if not os.path.exists(table_path):
        return "none"
    with open(table_path, newline="") as f:
        for row in csv.DictReader(f):
            if row["Model ID"] == model:
                return TYPE_TO_MODE.get(row["Type"], "none")
    return "none"


def top_similar(input_mol, cells: list[str], n: int) -> list[tuple[float, Chem.Mol]]:
    """The `n` distinct valid outputs most similar to `input_mol`, as
    (Tanimoto, mol) pairs, most similar first."""
    input_fp = FP_GEN.GetFingerprint(input_mol)
    scored = {}
    for cell in cells:
        if cell.strip().lower() in NULL_STRINGS:
            continue
        mol = Chem.MolFromSmiles(cell.strip())
        if mol is None:
            continue
        key = Chem.MolToSmiles(mol)
        if key not in scored:
            sim = DataStructs.TanimotoSimilarity(input_fp, FP_GEN.GetFingerprint(mol))
            scored[key] = (sim, mol)
    ranked = sorted(scored.items(), key=lambda item: (-item[1][0], item[0]))
    return [pair for _, pair in ranked[:n]]


def scaffold_highlights(input_mol, outputs: list) -> tuple[list[int], list[list[int]]]:
    """Atoms of the input's Murcko scaffold, and of each output that contains it."""
    scaffold = MurckoScaffold.GetScaffoldForMol(input_mol)
    if scaffold.GetNumAtoms() == 0:
        return [], [[] for _ in outputs]
    return list(input_mol.GetSubstructMatch(scaffold)), [
        list(mol.GetSubstructMatch(scaffold)) for mol in outputs
    ]


def ring_systems(mol) -> list[set[int]]:
    """Fused ring systems of `mol`, as sets of atom indices."""
    systems: list[set[int]] = []
    for ring in mol.GetRingInfo().AtomRings():
        merged = set(ring)
        for other in [s for s in systems if s & merged]:
            merged |= other
            systems.remove(other)
        systems.append(merged)
    return systems


def fragment_highlights(input_mol, outputs: list) -> tuple[list[int], list[list[int]]]:
    """Per output, the fragment it keeps from the input: the complete ring
    systems inside its largest common substructure with the input (the whole
    common substructure if there is none). For the input, the fragments kept
    by any of its outputs."""
    input_systems = ring_systems(input_mol)
    in_atoms: set[int] = set()
    out_atoms = []
    for mol in outputs:
        mcs = rdFMCS.FindMCS(
            [input_mol, mol],
            timeout=MCS_TIMEOUT_SECONDS,
            ringMatchesRingOnly=True,
            completeRingsOnly=True,
        )
        query = Chem.MolFromSmarts(mcs.smartsString) if mcs.numAtoms else None
        if query is None or mcs.numAtoms < MIN_FRAGMENT_ATOMS:
            out_atoms.append([])
            continue
        # atom q of the query is matched to in_match[q] and out_match[q]
        in_match = input_mol.GetSubstructMatch(query)
        out_match = mol.GetSubstructMatch(query)
        if not in_match or not out_match:
            out_atoms.append([])
            continue
        matched = set(in_match)
        kept = [
            q
            for q, atom in enumerate(in_match)
            if any(atom in system and system <= matched for system in input_systems)
        ] or list(range(len(in_match)))
        in_atoms |= {in_match[q] for q in kept}
        out_atoms.append([out_match[q] for q in kept])
    return sorted(in_atoms), out_atoms


def bonds_between(mol, atoms: list[int]) -> list[int]:
    """Bonds of `mol` whose two atoms are both in `atoms`."""
    chosen = set(atoms)
    return [
        b.GetIdx()
        for b in mol.GetBonds()
        if b.GetBeginAtomIdx() in chosen and b.GetEndAtomIdx() in chosen
    ]


def build_grid(inputs: list[str], output_rows: list[list[str]], mode: str) -> dict:
    """Everything the grid needs, cell by cell, row by row: molecules, legends,
    highlighted atoms and bonds, and the text to write into each blank cell
    (RDKit draws no legend for a missing molecule), keyed by cell position."""
    grid = {"mols": [], "legends": [], "atoms": [], "bonds": [], "blanks": {}}

    def add(mol, legend, atoms=()):
        grid["mols"].append(mol)
        grid["legends"].append(legend)
        grid["atoms"].append(list(atoms))
        grid["bonds"].append(bonds_between(mol, atoms) if mol is not None else [])

    for number, (smiles, cells) in enumerate(zip(inputs, output_rows), 1):
        input_mol = Chem.MolFromSmiles(smiles)
        if input_mol is None:
            grid["blanks"][len(grid["mols"])] = f"INPUT {number} (invalid SMILES)"
            for _ in range(N_SIMILAR + 1):
                add(None, "")
            continue
        best = top_similar(input_mol, cells, N_SIMILAR)
        outputs = [mol for _, mol in best]
        if mode == "scaffold":
            in_atoms, out_atoms = scaffold_highlights(input_mol, outputs)
        elif mode == "fragment":
            in_atoms, out_atoms = fragment_highlights(input_mol, outputs)
        else:
            in_atoms, out_atoms = [], [[] for _ in outputs]
        add(input_mol, f"INPUT {number}", in_atoms)
        for rank, ((sim, mol), atoms) in enumerate(zip(best, out_atoms), 1):
            add(mol, f"TS={sim:.2f} (top-{rank})", atoms)
        for _ in range(N_SIMILAR - len(best)):
            grid["blanks"][len(grid["mols"])] = "no output"
            add(None, "")
    return grid


def label_blanks(image, blanks: dict[int, str]) -> None:
    """Write the text of each blank cell where a legend would be drawn."""
    draw = ImageDraw.Draw(image)
    font = ImageFont.load_default(size=LEGEND_FONT_SIZE)
    for index, text in blanks.items():
        row, col = divmod(index, N_SIMILAR + 1)
        x = col * CELL_PIXELS + CELL_PIXELS // 2
        y = row * CELL_PIXELS + int(CELL_PIXELS * 0.92)
        draw.text((x, y), text, fill=(150, 150, 150), font=font, anchor="mm")


def table_row(model: str, table_path: str) -> dict:
    """The model's row of the master table, or an empty dict."""
    if os.path.exists(table_path):
        with open(table_path, newline="") as f:
            for row in csv.DictReader(f):
                if row["Model ID"] == model:
                    return row
    return {}


# Run with eos1klk's own interpreter: its t-SNE surrogate (Morgan fingerprint -> XGBoost).
TSNE_SNIPPET = """
import csv, sys
import numpy as np
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
from lazychemvis.artifacts.tsne import TSNEArtifact
checkpoints, src, dst = sys.argv[1:4]
smiles = [row[0] for row in list(csv.reader(open(src)))[1:]]
valid = [Chem.MolFromSmiles(smi) is not None for smi in smiles]
good = [smi for smi, ok in zip(smiles, valid) if ok]
coords = TSNEArtifact(dir_name=checkpoints).transform(good) if good else np.empty((0, 2))
rows = iter(coords)
with open(dst, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["tsne_x", "tsne_y"])
    for ok in valid:
        writer.writerow(next(rows) if ok else ["nan", "nan"])
"""


def project(smiles: list[str], models_dir: str, envs_dir: str) -> np.ndarray:
    """t-SNE coordinates (n x 2, NaN for unparseable SMILES) of `smiles` from
    eos1klk's t-SNE surrogate, run with that model's own env interpreter."""
    checkpoints = os.path.join(models_dir, PROJECTOR, "model", "checkpoints")
    python = os.path.join(envs_dir, f"{PROJECTOR}-cpu", "bin", "python")
    for label, path in (("checkpoints", checkpoints), ("env", python)):
        if not os.path.exists(path):
            raise SystemExit(
                f"{PROJECTOR} {label} not found: {path}; run 02_fetch_generative_models.py "
                f"--models {PROJECTOR} and 03_build_cpu_envs.py --models {PROJECTOR}"
            )
    env = {**os.environ, "PYTHONNOUSERSITE": "1", "PYTHONUNBUFFERED": "1"}
    with tempfile.TemporaryDirectory() as tmp:
        in_path, out_path = os.path.join(tmp, "in.csv"), os.path.join(tmp, "out.csv")
        with open(in_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["smiles"])
            writer.writerows([[smi] for smi in smiles])
        result = subprocess.run(
            [python, "-c", TSNE_SNIPPET, checkpoints, in_path, out_path],
            env=env,
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode != 0 or not os.path.exists(out_path):
            raise SystemExit(f"{PROJECTOR} failed: {result.stderr.strip()[-500:]}")
        with open(out_path, newline="") as f:
            rows = list(csv.reader(f))[1:]
    return np.array([[float(x), float(y)] for x, y in rows]).reshape(-1, 2)


def reference_coordinates(cache_dir: str, models_dir: str, envs_dir: str) -> np.ndarray:
    """t-SNE coordinates (m x 2) of the 20,000 random reference-library
    compounds, computed once through eos1klk and cached."""
    out_dir = os.path.join(cache_dir, PROJECTOR)
    cached = os.path.join(out_dir, f"reference_{N_REFERENCE // 1000}k.csv")
    if os.path.exists(cached):
        return np.loadtxt(
            cached, delimiter=",", skiprows=1, usecols=(1, 2), comments=None
        )
    os.makedirs(out_dir, exist_ok=True)
    library = os.path.join(cache_dir, "lazy-chemvis", "smiles_100k.csv")
    if not os.path.exists(library):
        os.makedirs(os.path.dirname(library), exist_ok=True)
        urllib.request.urlretrieve(REFERENCE_URL, library)
    smiles = [row[0] for row in read_rows(library) if row]
    sample = random.Random(REFERENCE_SEED).sample(smiles, N_REFERENCE)
    print(
        f"projecting {len(sample)} reference compounds with the {PROJECTOR} t-SNE model (once)"
    )
    coords = project(sample, models_dir, envs_dir)
    kept = [(smi, *xy) for smi, xy in zip(sample, coords) if not np.isnan(xy).any()]
    with open(cached, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["smiles", "tsne_x", "tsne_y"])
        writer.writerows(kept)
    print(f"  {len(kept)} of {len(sample)} reference compounds projected")
    return np.array([[x, y] for _, x, y in kept])


def seed_outputs(output_rows: list[list[str]]) -> list[list[str]]:
    """Per seed, up to MAX_OUTPUTS_PER_SEED distinct valid output SMILES."""
    per_seed = []
    for cells in output_rows:
        seen, kept = set(), []
        for cell in cells:
            smi = cell.strip()
            if smi.lower() in NULL_STRINGS:
                continue
            mol = Chem.MolFromSmiles(smi)
            if mol is None or Chem.MolToSmiles(mol) in seen:
                continue
            seen.add(Chem.MolToSmiles(mol))
            kept.append(smi)
        per_seed.append(kept[:MAX_OUTPUTS_PER_SEED])
    return per_seed


def point_density(xy: np.ndarray, low, high) -> np.ndarray:
    """Local density at each row of `xy`, scaled to 0-1: a 2D histogram over
    the plot extent, Gaussian-smoothed, read back at the point's cell."""
    counts, x_edges, y_edges = np.histogram2d(
        xy[:, 0],
        xy[:, 1],
        bins=DENSITY_BINS,
        range=[[low[0], high[0]], [low[1], high[1]]],
    )
    offsets = np.arange(-int(3 * DENSITY_SIGMA), int(3 * DENSITY_SIGMA) + 1)
    kernel = np.exp(-0.5 * (offsets / DENSITY_SIGMA) ** 2)
    kernel /= kernel.sum()
    smooth = np.apply_along_axis(np.convolve, 0, counts, kernel, "same")
    smooth = np.apply_along_axis(np.convolve, 1, smooth, kernel, "same")
    ix = np.clip(np.digitize(xy[:, 0], x_edges) - 1, 0, DENSITY_BINS - 1)
    iy = np.clip(np.digitize(xy[:, 1], y_edges) - 1, 0, DENSITY_BINS - 1)
    density = smooth[ix, iy]
    return density / density.max()


def dot_areas(xy: np.ndarray, sizes, low, high) -> np.ndarray:
    """Dot area of each row of `xy`, growing with its local density between
    the two `sizes`."""
    if len(xy) == 0:
        return np.empty(0)
    return sizes[0] + (sizes[1] - sizes[0]) * point_density(xy, low, high)


def draw_dots(ax, xy, areas, colour, alpha, zorder=1):
    """Scatter `xy` in one colour, the largest dots first so none hides another."""
    order = np.argsort(-areas)
    ax.scatter(
        xy[order, 0],
        xy[order, 1],
        s=areas[order],
        color=colour,
        alpha=alpha,
        linewidths=0,
        zorder=zorder,
    )


def plot_space(args, inputs: list[str], output_rows: list[list[str]]) -> str:
    """Draw the chemical-space map of the seeds and their outputs and return
    the path of the image."""
    reference = reference_coordinates(
        args.path_to_cache, args.path_to_models, args.path_to_envs
    )
    outputs = seed_outputs(output_rows)
    unique = list(dict.fromkeys(inputs + [smi for group in outputs for smi in group]))
    coords = dict(zip(unique, project(unique, args.path_to_models, args.path_to_envs)))

    def points(group: list[str]) -> np.ndarray:
        found = np.array([coords[smi] for smi in group]).reshape(-1, 2)
        return found[~np.isnan(found).any(axis=1)]

    overlay = [points(inputs)] + [points(group) for group in outputs]
    everything = np.vstack([reference] + [o for o in overlay if len(o)])
    margin = 0.03 * (everything.max(axis=0) - everything.min(axis=0))
    low, high = everything.min(axis=0) - margin, everything.max(axis=0) + margin

    os.environ.setdefault("MPLBACKEND", "Agg")  # stylia imports pyplot; run headless
    with mock.patch("shutil.rmtree"):  # stylia wipes matplotlib's cache dir on import
        import stylia  # only the map needs it

    # Format: print | Style: article; change with stylia.set_format() / set_style()
    stylia.set_format("print")
    stylia.set_style("article")
    palette = stylia.NamedColors()
    label_size = stylia.FONTSIZE_BIG
    reference_sizes = [f * stylia.MARKERSIZE for f in REFERENCE_SIZES]
    output_sizes = [f * stylia.MARKERSIZE for f in OUTPUT_SIZES]

    fig, axs = stylia.create_figure(1, len(inputs), height=MAP_HEIGHT)
    seed_xy = np.array([coords[smi] for smi in inputs])
    reference_areas = dot_areas(reference, reference_sizes, low, high)  # same in all
    for k, ((x, y), found) in enumerate(zip(seed_xy, overlay[1:])):
        ax = axs.next()
        draw_dots(ax, reference, reference_areas, palette.silver, 0.5)
        draw_dots(
            ax,
            found,
            dot_areas(found, output_sizes, low, high),
            palette.crimson,
            1,
            zorder=2,
        )
        if not np.isnan(x):
            ax.scatter(
                [x],
                [y],
                s=SEED_SIZE * stylia.MARKERSIZE,
                color=palette.black,
                linewidths=0,
                zorder=4,
            )
            ax.annotate(
                str(k + 1),
                (x, y),
                xytext=(3, 3),
                textcoords="offset points",
                fontsize=label_size,
                fontweight="bold",
                zorder=5,
            )
        ax.set_xlim(low[0], high[0])
        ax.set_ylim(low[1], high[1])
        ax.set_box_aspect(1)
        ax.grid(False)  # stylia switches the grid on
        ax.set_xticks([])
        ax.set_yticks([])
        stylia.label(ax, xlabel="tSNE-1", ylabel="tSNE-2" if k == 0 else "")
        ax.xaxis.label.set_fontsize(label_size)
        ax.yaxis.label.set_fontsize(label_size)
    row = table_row(args.model, args.table)
    title = " | ".join(x for x in (args.model, row.get("Slug"), row.get("Type")) if x)
    heading = fig.suptitle(
        title, fontsize=TITLE_FONT_FACTOR * label_size, verticalalignment="bottom"
    )
    fig.tight_layout()  # leaves room for the title
    first, last = fig.axes[0].get_position(), fig.axes[-1].get_position()
    heading.set_position(
        ((first.x0 + last.x1) / 2, first.y1 + TITLE_GAP_INCHES / fig.get_figheight())
    )
    out_dir = os.path.join(args.path_to_output, args.model, "png")
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, "chemical_space.png")
    stylia.save_figure(path)
    missing = sum(len(g) for g in outputs) + len(inputs) - sum(len(o) for o in overlay)
    print(
        f"{len(inputs)} seeds, {sum(len(o) for o in overlay[1:])} mapped outputs"
        f"{f' ({missing} could not be projected)' if missing else ''} -> {path}"
    )
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", required=True, help="Model identifier, e.g. eos2401")
    parser.add_argument(
        "--outputs",
        default=None,
        help="CSV of outputs for the example compounds (default: the model's "
        "shipped examples/run_output.csv)",
    )
    parser.add_argument(
        "--highlight",
        choices=["auto", "scaffold", "fragment", "none"],
        default="auto",
        help="What to highlight (default auto: from the model's Type in --table)",
    )
    parser.add_argument(
        "--plots",
        choices=["both", "grid", "space"],
        default="both",
        help="Which figure(s) to draw (default: both)",
    )
    parser.add_argument("--table", default=DEFAULT_TABLE)
    parser.add_argument("--path-to-models", default=DEFAULT_MODELS_DIR)
    parser.add_argument("--path-to-envs", default=DEFAULT_ENVS_DIR)
    parser.add_argument("--path-to-cache", default=DEFAULT_CACHE_DIR)
    parser.add_argument("--path-to-output", default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()

    examples = os.path.join(
        args.path_to_models, args.model, "model", "framework", "examples"
    )
    inputs_path = os.path.join(examples, "run_input.csv")
    outputs_path = args.outputs or os.path.join(examples, "run_output.csv")
    for label, path in (("inputs", inputs_path), ("outputs", outputs_path)):
        if not os.path.exists(path):
            raise SystemExit(f"{label} not found: {path}")

    inputs = [row[0] for row in read_rows(inputs_path) if row]
    output_rows = read_rows(outputs_path)
    if len(output_rows) != len(inputs):
        raise SystemExit(
            f"{outputs_path}: {len(output_rows)} rows for {len(inputs)} example compounds"
        )

    if args.plots in ("both", "space"):
        plot_space(args, inputs, output_rows)
    if args.plots == "space":
        return

    mode = args.highlight
    if mode == "auto":
        mode = mode_from_table(args.model, args.table)
    print(f"{args.model}: highlight = {mode}")

    grid = build_grid(inputs, output_rows, mode)
    options = Draw.MolDrawOptions()
    options.highlightBondWidthMultiplier = HIGHLIGHT_BOND_WIDTH
    options.highlightRadius = HIGHLIGHT_RADIUS
    options.legendFontSize = LEGEND_FONT_SIZE
    if mode in MODE_COLOUR:
        options.setHighlightColour(MODE_COLOUR[mode])
    image = Draw.MolsToGridImage(
        grid["mols"],
        molsPerRow=N_SIMILAR + 1,
        subImgSize=(CELL_PIXELS, CELL_PIXELS),
        legends=grid["legends"],
        highlightAtomLists=grid["atoms"],
        highlightBondLists=grid["bonds"],
        drawOptions=options,
    )
    label_blanks(image, grid["blanks"])

    out_dir = os.path.join(args.path_to_output, args.model, "png")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "examples_grid.png")
    image.save(out_path)

    per_row = N_SIMILAR + 1
    for i in range(len(inputs)):
        row = grid["legends"][i * per_row : (i + 1) * per_row]
        shown = [x for x in row[1:] if x]
        lit = sum(bool(a) for a in grid["atoms"][i * per_row + 1 : (i + 1) * per_row])
        best = f", best {shown[0]}" if shown else ""
        marked = f", highlight in {lit}" if mode != "none" else ""
        print(f"  compound {i + 1}: {len(shown)} of {N_SIMILAR} outputs{best}{marked}")
    print(f"{len(inputs)} x {per_row} grid -> {out_path}")


if __name__ == "__main__":
    main()
