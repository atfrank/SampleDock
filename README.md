# SampleDock
Molecular design framework the merges generative AI and molecular docking.

## Installation:
Anaconda is recommended to install SampleDock. However, the general requirement is
```
  - python=3.7
  - pytorch
  - cudatoolkit (optional only for cuda enabled device)
  - scipy
  - rdkit
  - rxdock (or locally compiled rdock)
```
To install the required python packages with anaconda:
```
conda env create -f environment.yml
source activate SampleDock
python setup.py install
```
If you are using Anaconda, the docking simulation program, *rDock*, is by default installed with the package `rxdock` (precompiled bin files). If you wish to compile it locally, it can be installed following the instruction: http://rdock.sourceforge.net/installation/. 

**Note**: the command options are slightly different between the two versions.

To run Sample and Dock, first specify the hyperparameters at `SnD/hyper.param`. Then, `python sampler.py`

## Software:

## Target Specific Libraries:
SARS-CoV-2: https://github.com/atfrank/SARS-CoV-2

