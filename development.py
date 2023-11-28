from chemsampler.samplers.sampler import UnitSampler

us = UnitSampler(model_id="eos1noy", timeout_sec=1)

smiles = us.get_example_smiles()
print(smiles)

print(us.sample(smiles))
