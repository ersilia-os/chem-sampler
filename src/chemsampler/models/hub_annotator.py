import pandas as pd

from ..hub.client import Backend, HubModel
from ..utils.logging import logger

#: A single oversized `run()` call degrades the served model and fails silently
#: past this scale (observed: 100k in one call collapsed to ~21% non-null after
#: ~25k molecules; six successive 10k calls against the same served model held
#: 99.3-99.5% non-null with zero errors). Chunk at this size regardless of input.
_MAX_CHUNK = 10_000


class HubAnnotator:
    """
    Scorer backed by an annotation model from the Ersilia Model Hub.

    Satisfies the same interface as `QEDAnnotator`, so it can be passed wherever an
    annotator is expected. Candidates are sent to the model in chunks of at most
    `_MAX_CHUNK`, each via its own `run()` call: these models are start-up dominated
    at the scale of a single round (hundreds to low thousands of molecules), so more
    calls cost little there, but a single call over tens of thousands of molecules
    has been observed to degrade the served model and fail silently.

    Parameters
    ----------
    model_id : str
        Ersilia identifier of the annotation model (e.g. "eos4zfy").
    column : str, optional
        Output column to score on. Defaults to the model's single numeric output
        column; required when the model returns more than one.
    backend : {"ersilia", "run_sh"}, optional
        Passed through to `HubModel`, by default "run_sh".
    """

    def __init__(
        self, model_id: str, column: str | None = None, backend: Backend = "run_sh"
    ):
        self.model_id = model_id
        self.column = column
        self._hub_model = HubModel(model_id, backend=backend)

    def score(self, smiles_list: list[str]) -> dict[str, float]:
        """
        Score a list of SMILES with the Hub model.

        Parameters
        ----------
        smiles_list : list[str]
            SMILES strings to score.

        Returns
        -------
        dict[str, float]
            Mapping from input SMILES to its score. Molecules the model could not
            process are omitted, so the result may be shorter than the input.
        """
        if not smiles_list:
            return {}

        scores = {}
        for start in range(0, len(smiles_list), _MAX_CHUNK):
            chunk = smiles_list[start : start + _MAX_CHUNK]
            df = self._hub_model.run(chunk)
            column = self.column or self._infer_column(df)
            for smi, value in zip(chunk, df[column]):
                if pd.notna(value):
                    scores[smi] = float(value)

        missing = len(smiles_list) - len(scores)
        if missing:
            logger.warning(
                f"{self.model_id}: {missing} of {len(smiles_list)} molecules scored null"
            )
        return scores

    def _infer_column(self, df: pd.DataFrame) -> str:
        """Pick the model's single numeric output column, ignoring input/key columns."""
        candidates = [
            c
            for c in df.columns
            if c not in ("key", "input", "smiles")
            and pd.api.types.is_numeric_dtype(df[c])
        ]
        if len(candidates) != 1:
            raise ValueError(
                f"{self.model_id} returns {len(candidates)} numeric columns {candidates}; "
                "pass column= to choose one."
            )
        return candidates[0]
