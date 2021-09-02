[![Anaconda-Server Badge](https://anaconda.org/aaronfranklab/sampledock/badges/installer/conda.svg)](https://conda.anaconda.org/aaronfranklab)
# SampleDock
* Molecular design framework that merges generative AI and molecular docking ([manuscript](https://www.biorxiv.org/content/10.1101/2020.06.09.143289v1) and [webpage](https://atfrank.github.io/SampleDock/))


## Installation:
### Option 1
Installation using Anaconda 

**Note:** The sampledock conda package is only distributed for linux systems. You will need to use **Option 2** for installation on Mac OS and Windows
```
conda create -n sampledock
conda activate sampledock
conda install sampledock -c AaronFrankLab -c conda-forge -c bioconda -c pytorch -c tmap
```
### Option 2
If you wish to install from source, the general requirement is listed below.
#### Main dependencies:
```
  - python=3.7
  - pytorch
  - cudatoolkit (optional only for cuda enabled device)
  - scipy
  - rdkit
  - rxdock (or rdock)
```
#### Data visualizations tools:
```
  - tmap
  - faerun
  - mhfp
```
You can clone this repo and install the required python packages with conda environment.yml
```
conda env create -f environment.yml
conda activate sampledock
python setup.py install
```
If you are using Anaconda and install the required packages with the `environment.yml`, the docking program, *rDock*, will be installed as the package `rxdock`. This is precompiled executables released by [RxDock](https://www.rxdock.org/). `rxdock` has some commandline argument changes from rDock. The original rDock can be installed following the instruction: http://rdock.sourceforge.net/installation/. 

## Usage:
To run Sample and Dock, first specify the hyperparameters in `hyper.param`. 

Then, `python -m sampledock hyper.param`

## JTVAE Model:
The JTVAE model is developed by Jin, W., Jaakkola, T. et al. and retireved from git: https://github.com/wengong-jin/icml18-jtnn. The default JTVAE generation model, `moses-h450z56`, is supplied with and trained with [MOSES dataset](https://github.com/molecularsets/moses).

The python scripts for JTVAE are modified for compatiblity with python 3.7.

## Target Specific Libraries:
SARS-CoV-2: https://github.com/atfrank/SARS-CoV-2


## COMMERCIAL USE LICENSE: 

If you are interested in commercial licensing of these applications (clinical, operational, etc.) please contact the University of Michigan Office of Technology Transfer for a quote and licensing options.

Drew Bennett - https://techtransfer.umich.edu/team/drew-bennett/

or

techtransfer@umich.edu
