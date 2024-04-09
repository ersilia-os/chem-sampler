from typing import List
from rdkit import Chem
from rdkit.Chem import DataStructs
import pandas as pd
import numpy as np

from ..samplers.sampler import UnitSampler
from ..descriptors.descriptor import DescriptorCalculator
from ..rules.rule import Ruler

class MasterSampler(object):
    def __init__(
        self,
        sampler_ids: List,
        descriptor_ids: List,
        unit_timeout_sec=60,
    ):
        self.sampler_ids = sampler_ids
        self.descriptor_ids = descriptor_ids
        self.unit_timeout_sec = unit_timeout_sec
        self._fetch()

    #automatically check if specified models are fetched, if not they will be fetched from Docker
    def _fetch(self):
        for sampler_id in self.sampler_ids:
            UnitSampler(model_id=sampler_id).fetch()
        for descriptor_id in self.descriptor_ids:
            DescriptorCalculator(model_id=descriptor_id).fetch()
    
    def _clean_sampled_smiles(self, sampled_smiles, keep_smiles=None, avoid_smiles=None):
        rule = Ruler(keep_smiles, avoid_smiles)
        sampled_smiles_filtered = sampled_smiles
        if keep_smiles is not None:
            sampled_smiles_filtered = rule.keep_substructure(sampled_smiles)
        print("KEEP SUBSTRUCTURE", len(sampled_smiles_filtered))
        if len(sampled_smiles_filtered)==0:
            print("No SMILES keeps the core structure")
            return None
        if avoid_smiles is not None:
            sampled_smiles_filtered = rule.avoid_substructure(sampled_smiles_filtered)
        print("AVOID SUBSTRUCTURE", len(sampled_smiles_filtered))
        return sampled_smiles_filtered
    
    def _sample(self, input_smiles, keep_smiles=None, avoid_smiles=None):
        sampled_smiles_all = set()
        sampler_info = {}
        for sampler_id in self.sampler_ids:
            us = UnitSampler(model_id=sampler_id, timeout_sec=self.unit_timeout_sec)
            sampled_smiles = us.sample(input_smiles)
            print(sampler_id, len(sampled_smiles))
            sampler_info[sampler_id]=[len(sampled_smiles)]
            if keep_smiles is not None or avoid_smiles is not None:
                sampled_smiles = self._clean_sampled_smiles(sampled_smiles, keep_smiles, avoid_smiles)
            if sampled_smiles is not None:
                sampler_info[sampler_id].append(len(sampled_smiles))
                sampled_smiles_all.update(sampled_smiles)
            else: 
                sampler_info[sampler_id].append(0)
            print("TOTAL SMILES", len(sampled_smiles_all))
        sampled_smiles_all = list(sampled_smiles_all)
        return sampled_smiles_all, sampler_info
    
    def _calculate_seed_descriptors(self, seed_smiles, descriptor_id):
        dc = DescriptorCalculator(model_id=descriptor_id)
        descs = dc.calculate([seed_smiles])[0]
        return descs
        
    def _calculate_sampled_descriptors(self, sampled_smiles, descriptor_id):
        dc = DescriptorCalculator(model_id=descriptor_id)
        print(sampled_smiles)
        descs = dc.calculate(sampled_smiles)
        return descs

    def _calculate_euclidean_distance (self, a1, a2):
        euds = np.linalg.norm(a1 - a2)
        return euds
    
    def _calculate_tanimoto_similarity(self, fp1, fp2):
        return DataStructs.TanimotoSimilarity(fp1, fp2)

    def _np_to_bv(self, fv):
        bv = DataStructs.ExplicitBitVect(len(fv))
        for i,v in enumerate(fv):
            if v:
                bv.SetBit(i)
        return bv

    def _check_descriptor_output_type(self,descriptor_id):
        dc = DescriptorCalculator(descriptor_id)
        info = dc.get_info()
        return info["metadata"]["Output Type"][0]

    def _check_descriptor_sparse(self,desc):
        zeroes = np.count_nonzero(desc == 0) / len(desc)
        if zeroes > 0.5:
            return True
        else:
            return False

    def _is_binary(self,desc):
        unique_values = np.unique(desc)
        return np.array_equal(unique_values, np.array([0, 1]))

    def _calculate_similarities(self, seed_smiles, sampled_smiles):
        similarities_dict = {}
        for descriptor_id in self.descriptor_ids:
            origin_descs = self._calculate_seed_descriptors(seed_smiles, descriptor_id)
            sampled_descs = self._calculate_sampled_descriptors(sampled_smiles, descriptor_id)
            if self._check_descriptor_output_type(descriptor_id) == "Float":
                similarities = [self._calculate_euclidean_distance(origin_descs, sampled_desc) for sampled_desc in sampled_descs]
                similarities_dict[descriptor_id+"_euclidean"] = similarities              
            elif self._check_descriptor_output_type(descriptor_id) == "Integer":
                if self._check_descriptor_sparse(origin_descs) is False:
                    similarities = [self._calculate_euclidean_distance(origin_descs, sampled_desc) for sampled_desc in sampled_descs]
                    similarities_dict[descriptor_id+"_euclidean"] = similarities  
                else:
                    if self._is_binary(origin_descs) is True:
                        similarities = [self._calculate_tanimoto_similarity(self._np_to_bv(origin_descs), self._np_to_bv(sampled_desc)) for sampled_desc in sampled_descs]
                        similarities_dict[descriptor_id+"_tanimoto"] = similarities  
                    else:
                        similarities = [self._calculate_euclidean_distance(origin_descs, sampled_desc) for sampled_desc in sampled_descs]
                        similarities_dict[descriptor_id+"_euclidean"] = similarities                          
            else:
                print("Output Type unknown")
        return similarities_dict

    def run(self, seed_smiles, input_smiles=None, keep_smiles=None, avoid_smiles=None):
        if input_smiles == None:
            input_smiles = seed_smiles #in the first round, seed and input are the same
        sampled_smiles, sampled_info = self._sample(input_smiles, keep_smiles, avoid_smiles)
        print("SAMPLED", len(sampled_smiles))
        sampled_smiles = [smi for smi in sampled_smiles if smi is not None]
        print("ALL VALID SAMPLED", len(sampled_smiles))
        similarities_dict = self._calculate_similarities(seed_smiles, sampled_smiles)
        df = pd.DataFrame()
        df["sampled_smiles"] = sampled_smiles
        for k,v in similarities_dict.items():
            df[k] = v
        return df, sampled_info
    

