from chemsampler.samplers.bimodal.sampler import BimodalSampler
n =  10
smiles_list = []
samples = BimodalSampler().sample(smiles_list ,n)
samples.to_csv('./bimodal_smiles.csv', index=False)
