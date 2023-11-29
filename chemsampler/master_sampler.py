import random
from .samplers.sampler import UnitSampler
from .descriptors.descriptor import DescriptorCalculator


class MasterSampler(object):
    def __init__(
        self,
        sampler_ids,
        descriptor_ids,
        num_samples,
        num_rounds=100,
        unit_timeout_sec=60,
    ):
        self.sampler_ids = sampler_ids
        self.descriptor_ids = descriptor_ids
        self.num_samples = num_samples
        self.unit_timeout_sec = unit_timeout_sec
        self.num_rounds = num_rounds
        self._fetch()

    def _fetch(self):
        for sampler_id in self.sampler_ids:
            UnitSampler(model_id=sampler_id).fetch()
        for descriptor_id in self.descriptor_ids:
            DescriptorCalculator(model_id=descriptor_id).fetch()

    def _precalculate_input_descriptors(self, smiles):
        data = {}
        for descriptor_id in self.descriptor_ids:
            data[descriptor_id] = DescriptorCalculator(
                model_id=descriptor_id
            ).calculate([smiles])[0]
        self.input_descriptors = data

    def run(self, smiles):
        self._precalculate_input_descriptors(smiles)
        sampled_smiles = set()
        for sampler_id in self.sampler_ids:
            us = UnitSampler(model_id=sampler_id, timeout_sec=self.unit_timeout_sec)
            sampled_smiles.update(us.sample(smiles))
        sampled_smiles = list(sampled_smiles)
        for descriptor_id in self.descriptor_ids:
            dc = DescriptorCalculator(model_id=descriptor_id)
            descs = dc.calculate(sampled_smiles)
