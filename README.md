# ChemSampler
Sample small molecules, both from large chemical libraries as well as generative models

## Installation

It requires python>=3.7. Create a conda environment and activate it
```python
conda create -n chemsampler pthon=3.7
conda activate chemsampler
```
Clone the github repo:

```
git clone https://github.com/ersilia-os/chem-sampler.git
cd chem-sampler
pip install -e .
```

# chemsampler class 
It samples from large chemical libraries and pre-sampled molecules in chemsample database from generative models.

## Usage

```python
from chemsampler import ChemSampler
from chemsampler import example

smiles_list = example()
sampler = ChemSampler()
sampled_smiles = sampler.sample(smiles_list, num_samples=1000, sim_ub=0.95, sim_lb=0.6, distribution="ramp")
```

# Generative Samplers

## 1. Moler 

Moler is a generative model of molecular graphs developed by microsoft research  which supports scaffold-constrained generation. Publication could be found at [Learning to Extend Molecular Scaffolds with Structural Motifs](https://arxiv.org/abs/2103.03864)

To use moler build conda environment :

```
cd chemsampler/tools/moler
conda env create -f environment.yml
```

### API available for Moler
### To Sample : sample
Sampling from pre-calculated molecules
```python
from chemsampler.samplers.moler.sampler import MolerSampler
smiles_list = ['C1(C2=CC=C(NC(=O)C(N)C)C=C2)=CSC3=C1C=CC=C3']
sampler = MolerSampler()
samples = MolerSampler().sample(n= 10, smiles=smiles_list, search_pre_calculated=True, cutoff= 0.5)
print(samples)

```
sampling from pre-trained model

```python
from chemsampler.samplers.moler.sampler import MolerSampler

sampler = MolerSampler()
samples = sampler.sample(n= 10)
print(samples)
```
### To fine-tune : fit 
fit api fine-tunes the data on custom dataset.Input to fit is as list of smiles and the file name to save the newly generated checkpoint. The same file name can be used to generate from the fitted model.

```python
from chemsampler.samplers.moler.fit import MolerFit
from chemsampler.samplers.moler.sampler import MolerSampler

from chemsampler import example
checkpoint_dir = 'fit_example_checkpoint_dir'
fit = MolerFit()
smiles_list = example()
fit(smiles_list, checkpoint_dir)
#sample from finetuned model
sampler = MolerSampler()
samples = sampler.sample(n= 100, checkpoint_dir= checkpoint_dir)
print(samples)
```

## 2. Bimodal 

Bimodal is Bidirectional Molecule Generation model with Recurrent Neural Networks developed by ETH modal lab. The publication can be found at [Learning to Extend Molecular Scaffolds with Structural Motifs](https://arxiv.org/abs/2103.03864) .

[Github](https://github.com/ETHmodlab/BIMODAL)

To use bimodal build conda environment :

```
cd chemsampler/tools/bimodal
conda env create -f environment.yml
```

### API available : sample anf fit
To sample from bimodal `from chemsampler.samplers.bimodal.sampler import BimodalSampler`. The sample works same as Moler.

To fit bimodal import `from chemsampler.samplers.bimodal.fit import BimodalFit`. The fit method works similar to Moler. 
It takes list of smiles as input for fint-tuning.

## 3. mgm(Masked Graph modeling)
This model has been developed by NYU Deel Learning lab. It's faster than the other generative models available. Publication can be found at [Masked graph modeling for molecule generation](https://www.nature.com/articles/s41467-021-23415-2)

Only sample API is availble for mgm. 
To use mgm build conda environment :

```
cd chemsampler/tools/mgm
conda env create -f environment.yml
```
then import mgm `from chemsampler.samplers.mgm.sampler import MgmSampler`

## 4. hybridCLM 

HybridCLM leverages molecular structure and bioactivity with chemical language models to sample molecules. Ideally this model should be fine tuned to focus on targeted bioactivity. M

[Publication](https://chemrxiv.org/engage/chemrxiv/article-details/615580ced1fc334326f9356e) and 
[Github](https://github.com/ETHmodlab/hybridCLMs)

To use this model build conda environment :

```
cd chemsampler/tools/hybridclm
conda env create -f environment.yml
```
### Fine-tune and sample
To fine tune this model smiles along with bioactivity is required.

## 5. fastjtnn
To be replaced by hgraph2grapg which is faster
