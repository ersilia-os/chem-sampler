from chemsampler.master_sampler import MasterSampler


ms = MasterSampler(sampler_ids=["eos1noy", "eos1d7r"], descriptor_ids = ["eos5axz", "eos1ut3"], unit_timeout_sec = 60)
seed_smiles = "FC(F)OC(C=C1)=CC=C1C2=NN=C3C=NC=C(OCC(CO)C4=C(F)C=CC=C4)N32"
keep_smarts = "[*]-c1cncc2nnc(-[*]-3-[#6]-[#6]-[*]-[#6]-[#6]-3)n12"
avoid_smarts = "F[O](F)[#6]-[#6]-1=[#6]-[#6]=[#6]-[#6]=[#6]-1"
df = ms.run(seed_smiles=seed_smiles, keep_smarts=None, avoid_smarts=None)
print(df.shape)
df.to_csv("out.csv")