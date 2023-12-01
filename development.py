from chemsampler.master_sampler import MasterSampler


ms = MasterSampler(sampler_ids=["eos1noy", "eos1d7r"], descriptor_ids = ["eos5axz", "eos1ut3"], unit_timeout_sec = 60)
input_smiles = "CN1CCN(C2=c3ccccc3=Nc3ccc(Cl)cc3N2)CC1"
df = ms.run(input_smiles)
print(df.head)