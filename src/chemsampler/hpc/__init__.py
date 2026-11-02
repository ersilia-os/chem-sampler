"""Running a whole `chemsampler run` on a Slurm cluster (`run --hpc NAME`).

Each round, the similarity and generative models run as parallel jobs; once all
are done their candidates are put together, then every annotator (Hub models and
RDKit ones) runs as its own job, and the round is ranked as a local run would.
The profiles are in `hpc.json` in chemsampler's folder (see `profiles`).
"""

from .executor import SlurmExecutor, prepare
from .profiles import Profile, load_profile

__all__ = ["Profile", "SlurmExecutor", "load_profile", "prepare"]
