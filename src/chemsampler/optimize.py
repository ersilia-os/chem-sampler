import pandas as pd

from .utils.logging import logger


def hill_climb(
    seed_smiles: str, generator, annotator, n_rounds: int = 5
) -> pd.DataFrame:
    """
    Iteratively improve a seed molecule's score by alternating generation and scoring.

    Each round, `generator` proposes candidates from the current-best molecule and
    `annotator` scores them. The best-scoring candidate becomes the seed for the
    next round. Stops after `n_rounds` or as soon as a round fails to beat the
    current best.

    Parameters
    ----------
    seed_smiles : str
        SMILES of the starting molecule.
    generator : object
        Any object exposing `.generate(smiles: str) -> list[str]`.
    annotator : object
        Any object exposing `.score(smiles_list: list[str]) -> dict[str, float]`.
    n_rounds : int, optional
        Maximum number of rounds to run, by default 5.

    Returns
    -------
    pandas.DataFrame
        One row per round, with columns `round`, `smiles`, `score`, `is_new_best`.
    """
    best_smiles = seed_smiles
    best_score = annotator.score([seed_smiles])[seed_smiles]
    history = [
        {"round": 0, "smiles": best_smiles, "score": best_score, "is_new_best": True}
    ]
    logger.info(f"Round 0: seed score = {best_score:.3f}")

    for round_num in range(1, n_rounds + 1):
        candidates = generator.generate(best_smiles)
        if not candidates:
            logger.warning(f"Round {round_num}: no candidates generated, stopping.")
            break

        scores = annotator.score(candidates)
        round_best_smiles = max(scores, key=scores.get)
        round_best_score = scores[round_best_smiles]
        is_new_best = round_best_score > best_score

        history.append(
            {
                "round": round_num,
                "smiles": round_best_smiles,
                "score": round_best_score,
                "is_new_best": is_new_best,
            }
        )

        if is_new_best:
            logger.success(f"Round {round_num}: improved to {round_best_score:.3f}")
            best_smiles, best_score = round_best_smiles, round_best_score
        else:
            logger.info(f"Round {round_num}: no improvement, stopping.")
            break

    return pd.DataFrame(history)
