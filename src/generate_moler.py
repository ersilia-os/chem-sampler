from chemsampler.samplers.moler.sampler import MolerSampler
sampler = MolerSampler()
n = 1000000
smiles_list = []
samples = MolerSampler().sample(smiles_list, n)
samples.to_csv('./moler_smiles.csv', index=False)

