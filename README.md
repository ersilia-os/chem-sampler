# ChemSampler
ChemSampler is a simple Python API to sample the chemical space around a given molecule. It is developed and maintained by the [Ersilia Open Source Initiative](https://ersilia.io).

ChemSampler capitalizes on the *generative* models available in the [Ersilia Model Hub](https://ersilia.io/model-hub). Feel free to suggest new models by [opening an issue](https://github.com/ersilia-os/ersilia/issues).

> A legacy version of ChemSampler is available in the archived [`chem-sampler-legacy`](https://github.com/ersilia-os/chem-sampler-legacy) repository.

## Warning
This repository is work-in-progress. We are actively working on it and hence new features and improved documentation will be updated on a weekly basis.

## Installation

Create a Conda environment and activate it:
```bash
conda create -n chemsampler python=3.10
conda activate chemsampler
```

Chem-Sampler only requires Ersilia and RdKit to run. Please make sure to install the latest version of the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia).

## Usage
ChemSampler allows the user to specify a list of samplers and a list of descriptors from the Ersilia Model Hub. Providing an initial input molecule, ChemSampler will generate new candidates based on the list of samplers and select the best molecules using similarity searches with the provided descriptors. For good results, users should evaulate the type of descriptors used, for example those related to pharmacophore and 3D structure, in addition to traditional Morgan and other 2D descriptors

## WIP
We are adding new features to improve the scaffold hopping capabilities of ChemSampler, particularly thinking about natural products. We will allow for a number of rules (i.e, do not maintain certain core structure) to be specified by the user. Check GitHub Projects to follow the development of new features

## License
All the code in this repository is licensed under a GPL-v3 License. Please note that individual models fetched through ChemSampler might have different open source licenses.

## About
Ersilia is a non-profit organisation developing AI tools for infectious disease research in the Global South. Learn [more](https://ersilia.io) about us and [Support our work!](https://ersilia.io/donate)

