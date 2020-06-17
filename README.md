# SampleDock
Molecular design framework the merges generative AI and molecular docking. https://www.biorxiv.org/content/10.1101/2020.06.09.143289v1

## Installation:
Anaconda is recommended to install SampleDock. However, the general requirement is listed below.
```
  - python=3.7
  - pytorch
  - cudatoolkit (optional only for cuda enabled device)
  - scipy
  - rdkit
  - rxdock (or rdock)
```
To install the required python packages with anaconda:
```
conda env create -f environment.yml
source activate sampledock
python setup.py install
```
If you are using Anaconda, the docking simulation program, *rDock*, is by default installed with the package `rxdock` (precompiled bin files). If you wish to compile it locally, it can be installed following the instruction: http://rdock.sourceforge.net/installation/. 

**Note**: the command options are slightly different between the two versions.

To run Sample and Dock, first specify the hyperparameters at `hyper.param`. 

Then, `python sampler.py -p hyper.param`

## JTVAE Model:
The JTVAE model is developed by Jin, W., Jaakkola, T. et al. and retireved from git: https://github.com/wengong-jin/icml18-jtnn. The default JTVAE generation model, `moses-h450z56`, is supplied with and trained with [MOSES dataset](https://github.com/molecularsets/moses).

The python script for JTVAE is updated for compatiblity with python 3.7.

## Target Specific Libraries:
SARS-CoV-2: https://github.com/atfrank/SARS-CoV-2

