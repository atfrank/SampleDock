# SampleDock
Molecular design framework the merges generative AI and molecular docking

## Installation:
Anaconda is recommended. To install the required python packages:
```
conda env create -f environment.yml
source activate SampleDock
python setup.py install
```
The docking simulation program, `rdock` is installed by conda with package rxdock (precompiled bin files). If you wish to compile it locally, it can be installed following the instruction: http://rdock.sourceforge.net/installation/. 

**Note**: command options are slightly different between the two versions.

To run Sample and Dock, first specify the hyperparameters at `SnD/hyper.param`. Then, `python sampler.py`

## Software:

## Target Specific Libraries:
SARS-CoV-2 https://github.com/atfrank/SARS-CoV-2

