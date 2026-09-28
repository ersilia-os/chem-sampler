"""
Fetch every generative model from the Ersilia Model Hub catalog.

Pulls the public catalog (https://catalog.ersilia.io/api/models), keeps
models where Task == "Sampling" and Subtask == "Generation", reports how
many there are, and for each one: clones its GitHub repo (or pulls it if
already cloned) into --path-to-models, then runs `eosvc download --path .`
inside it to fetch the model's checkpoint(s) from the public eosvc bucket.
--models overrides this: skips the catalog and processes exactly the given
identifier(s) instead.

Fully self-contained: no imports from chemsampler. Requires `git` and the
`eosvc` CLI (installed via the eosvc package) on PATH.
"""

import argparse
import json
import os
import subprocess
import urllib.request

CATALOG_API_URL = "https://catalog.ersilia.io/api/models"
GITHUB_ORG = "ersilia-os"

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_MODELS_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "models"))


def fetch_catalog() -> list[dict]:
    """Download the full Ersilia Model Hub catalog as a list of model records."""
    print(f"Fetching {CATALOG_API_URL}")
    with urllib.request.urlopen(CATALOG_API_URL) as resp:
        return json.load(resp)


def filter_generative(models: list[dict]) -> list[dict]:
    """Keep only Task == "Sampling" and Subtask == "Generation" models."""
    return [
        m
        for m in models
        if m.get("Task") == "Sampling" and m.get("Subtask") == "Generation"
    ]


def clone_or_pull(identifier: str, models_dir: str) -> bool:
    """Clone `identifier`'s repo into models_dir/, or pull it if already cloned.

    Returns True on success. The clone directory name must equal `identifier`
    exactly - eosvc infers the S3 path to download from from the current
    working directory's basename.
    """
    model_dir = os.path.join(models_dir, identifier)
    if os.path.isdir(os.path.join(model_dir, ".git")):
        print(f"  {identifier}: already cloned, pulling")
        result = subprocess.run(
            ["git", "-C", model_dir, "pull", "--quiet"],
            capture_output=True,
            text=True,
            check=False,
        )
    else:
        print(f"  {identifier}: cloning")
        repo_url = f"https://github.com/{GITHUB_ORG}/{identifier}.git"
        result = subprocess.run(
            ["git", "clone", "--quiet", repo_url, model_dir],
            capture_output=True,
            text=True,
            check=False,
        )
    if result.returncode != 0:
        print(f"  {identifier}: git failed - {result.stderr.strip()}")
        return False
    return True


def download_checkpoints(identifier: str, models_dir: str) -> bool:
    """Run `eosvc download --path .` from inside the model's own directory."""
    model_dir = os.path.join(models_dir, identifier)
    result = subprocess.run(
        ["eosvc", "download", "--path", "."],
        cwd=model_dir,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.stdout.strip():
        print(result.stdout.strip())
    if result.returncode != 0:
        print(f"  {identifier}: eosvc reported an issue (see above), continuing")
    return result.returncode == 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--path-to-models",
        default=DEFAULT_MODELS_DIR,
        help=f"Where to clone models into (default: {DEFAULT_MODELS_DIR})",
    )
    parser.add_argument(
        "--models",
        nargs="+",
        default=None,
        help="One or more model identifiers (e.g. eos9taz). Supersedes the "
        "catalog: skips fetching/filtering and downloads exactly these.",
    )
    args = parser.parse_args()
    os.makedirs(args.path_to_models, exist_ok=True)

    if args.models:
        generative = [{"Identifier": identifier} for identifier in args.models]
        print(f"Using {len(generative)} model(s) from --models, skipping catalog")
    else:
        generative = filter_generative(fetch_catalog())
        print(f"Found {len(generative)} generative (Sampling / Generation) models")

    n_cloned = 0
    n_checkpoints_ok = 0
    for model in generative:
        identifier = model["Identifier"]
        if clone_or_pull(identifier, args.path_to_models):
            n_cloned += 1
            if download_checkpoints(identifier, args.path_to_models):
                n_checkpoints_ok += 1

    print(
        f"Done: {n_cloned}/{len(generative)} repos ready, "
        f"{n_checkpoints_ok}/{len(generative)} checkpoint downloads confirmed"
    )


if __name__ == "__main__":
    main()
