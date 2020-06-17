import argparse

parser = argparse.ArgumentParser(prog='Sample-n-Dock')
parser.add_argument("params",metavar = 'P', type=str, 
                    help="parameter file for Sample-n-Dock")
parser.add_argument("-o","--output", help="directory for output, default to ../SnD_designs", 
                    default = "../SnD_designs")
a = parser.parse_args()

if a.params is None:
    print('A parameter file is needed! \n'+
          'Example can be found in: https://github.com/atfrank/SampleDock/blob/master/hyper.param')
    exit()
    
import torch
import os
import subprocess

from rdkit import rdBase
## Disable rdkit Logs
rdBase.DisableLog('rdApp.error')
from .jtvae import Vocab, JTNNVAE
# Sample and Dock tools
from .SnD import prep_prm
from .SnD import dock, sort_pose, save_pose
from .SnD import hyperparam_loader, create_wd, smiles_to_sdfile

# Load hyper parameters
p = hyperparam_loader(a.params)

## Load Stock JTNN VAE Model
vocab = Vocab(p.vocab)
jtvae = JTNNVAE(vocab, p.hidden_size, p.latent_size, p.depthT, p.depthG)

if torch.cuda.is_available():
    device = torch.device("cuda")
    print('Using CUDA device:',torch.cuda.get_device_name(torch.cuda.current_device()))
else:
    device = torch.device("cpu")
    print("Using CPU for torch device")

jtvae.load_state_dict(torch.load(p.model_loc, map_location=device))

## VAE encoding and decoding on the initial seeding smiles
try:
    design_list = jtvae.smiles_gen(p.seed_smi, p.ndesign)
except KeyError as err:
    print('[KeyError]',err,
          'does not exist in the current JTVAE model vocabulary "%s" (the training set of the model did not contain this structure),'%p.vocab_loc,
          'thus "%s" failed to initialize the model as seeding molecule!'%p.seed_smi)
    exit()

# Create working directory
wd = create_wd(a.output,p.receptor_name)

## prepare rdock .prm file and get file path
prmfile, cav_dir = prep_prm(p.receptor_file,p.ligand_file,p.receptor_name,wd)

## create pocket
cmdline = p.cavity_protocol+' %s > %s/create_cavity.out'%(prmfile,cav_dir)
proc = subprocess.Popen(cmdline, shell=True)
proc.wait()
print('Docking pocket grid created')

## Main loop: VAE on subsequent returned compounds
for j in range(p.ncycle):

    design_dir = os.path.abspath(os.path.join(wd,'cycle_%s'%j))
    docking_dir = os.path.abspath(os.path.join(design_dir, 'docking'))
    os.makedirs(docking_dir)

    ## write .sdf file and get ligs file names
    ligs = smiles_to_sdfile(design_list,design_dir)

    ## Generate docking scores from .sd files
    dock(ligs, docking_dir, prmfile, p.docking_prm, p.npose, p.prefix)
    ranked_poses = sort_pose(docking_dir, p.sort_by, p.prefix) 
    save_pose(ranked_poses, design_dir)

    ## Generate new design list
    for energy, name, mol in ranked_poses:
        smi = mol.GetProp('SMILES')
        design_list = []
        try:
            print('[INFO]: Generating new designs \t', end = '\r')
            design_list = jtvae.smiles_gen(smi, p.ndesign)
        # go to the second best candidate if the best does not give any return
        except KeyError as err:
            print('[KeyError]',err,'is not part of the vocabulary')
            continue

        if len(design_list) != 0: 
            break 

        else: 
            print('Current design (%s) has no offspring; trying the next one \r'%name)

    print("[INFO]: Cycle %s: %s %s kcal/mol"%(j, smi, energy)+'\t'*6)

