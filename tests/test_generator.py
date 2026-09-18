from chemsampler.models.generator import GeneratorPool


class StubGenerator:
    """Stands in for a HubGenerator without touching the Ersilia Model Hub."""

    def __init__(self, model_id, candidates):
        self.model_id = model_id
        self._candidates = candidates

    def generate(self, seed_smiles):
        return list(self._candidates)


def test_pool_unions_candidates_and_deduplicates():
    pool = GeneratorPool(
        [
            StubGenerator("model_a", ["CCO", "CCC"]),
            StubGenerator("model_b", ["CCC", "CCN"]),
        ]
    )

    candidates = pool.generate("CC")

    assert sorted(candidates) == ["CCC", "CCN", "CCO"]


def test_pool_keeps_provenance():
    pool = GeneratorPool(
        [
            StubGenerator("model_a", ["CCO", "CCC"]),
            StubGenerator("model_b", ["CCC"]),
        ]
    )

    by_model = pool.generate_by_model("CC")

    assert by_model == {"model_a": ["CCO", "CCC"], "model_b": ["CCC"]}
