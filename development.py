from chemsampler.descriptors.descriptor import DescriptorCalculator

dc = DescriptorCalculator(model_id="eos1ut3")
dc.fetch()

smiles_list = dc.get_example_smiles_list()
print(smiles_list)

print(dc.calculate(smiles_list))
