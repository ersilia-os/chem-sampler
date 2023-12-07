from chemsampler.master_sampler import MasterSampler


ms = MasterSampler(sampler_ids=["eos1noy", "eos8fma"], descriptor_ids = ["eos5axz", "eos1ut3"], unit_timeout_sec = 60)
seed_smiles = "FC(F)OC(C=C1)=CC=C1C2=NN=C3C=NC=C(OCC(CO)C4=C(F)C=CC=C4)N32"
keep_smiles = "*c1nnc2cncc(*)n12"
avoid_smiles = "F[O](F)CC1=CC=CC=C1"
df = ms.run(seed_smiles=seed_smiles, keep_smiles=None, avoid_smiles=None)
df.to_csv("~/out.csv", index=False)