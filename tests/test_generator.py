from chemsampler.models.generator import GeneratorPool, SeedRequiredError


class StubGenerator:
    """Stands in for a HubGenerator without touching the Ersilia Model Hub."""

    def __init__(self, model_id, candidates):
        self.model_id = model_id
        self._candidates = candidates

    def generate(self, seed_smiles):
        return list(self._candidates)


class SeedRequiredStubGenerator:
    """Stands in for a Hub generator that cannot run without a seed."""

    model_id = "seed_required"

    def generate(self, seed_smiles):
        if seed_smiles is None:
            raise SeedRequiredError(f"{self.model_id} requires a seed molecule")
        return ["CCO"]


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


def test_pool_warns_on_empty_generator(caplog):
    pool = GeneratorPool(
        [
            StubGenerator("model_a", ["CCO"]),
            StubGenerator("model_b", []),
        ]
    )

    with caplog.at_level("WARNING"):
        candidates = pool.generate("CC")

    assert candidates == ["CCO"]
    assert any(
        "model_b" in r.message and "0 candidates" in r.message for r in caplog.records
    )


def test_pool_by_model_warns_on_empty_generator(caplog):
    pool = GeneratorPool(
        [
            StubGenerator("model_a", ["CCO"]),
            StubGenerator("model_b", []),
        ]
    )

    with caplog.at_level("WARNING"):
        by_model = pool.generate_by_model("CC")

    assert by_model == {"model_a": ["CCO"], "model_b": []}
    assert any(
        "model_b" in r.message and "0 candidates" in r.message for r in caplog.records
    )


def test_pool_skips_generator_that_requires_a_seed(caplog):
    pool = GeneratorPool(
        [SeedRequiredStubGenerator(), StubGenerator("model_b", ["CCN"])]
    )

    with caplog.at_level("WARNING"):
        candidates = pool.generate(None)

    assert candidates == ["CCN"]
    assert any(
        "seed_required" in r.message and "requires a seed" in r.message
        for r in caplog.records
    )


def test_pool_by_model_skips_generator_that_requires_a_seed(caplog):
    pool = GeneratorPool(
        [SeedRequiredStubGenerator(), StubGenerator("model_b", ["CCN"])]
    )

    with caplog.at_level("WARNING"):
        by_model = pool.generate_by_model(None)

    assert by_model == {"seed_required": [], "model_b": ["CCN"]}
