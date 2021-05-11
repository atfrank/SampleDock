## TODO: this will definitely need to be restructured as the additional features have cluterred code.
import argparse

parser = argparse.ArgumentParser(prog='sampledock')
parser.add_argument("params",metavar = 'P', type=str, 
                    help="parameter file for SampleDock. For a parameter file template, please see \
                        https://raw.githubusercontent.com/atfrank/SampleDock/master/hyper.param")
parser.add_argument("-o","--output", help="directory for output, default to ../SnD_designs", 
                    default = "../SnD_designs")
a = parser.parse_args()

import torch
import os
import sys
import subprocess
import pickle
from rdkit import rdBase, Chem
## Disable rdkit Logs
rdBase.DisableLog('rdApp.error')
from .jtvae import Vocab, JTNNVAE
# Sample and Dock tools
from .SnD import prep_prm
from .SnD import dock, sort_pose, save_pose
from .SnD import hyperparam_loader, create_wd, smiles_to_sdfile
from .SnD import single_generator, distributed_generator
from .SnD import combine_designs, mkdf
from .SnD import LSH_Convert, tree_coords, df_to_faerun
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
    print("CUDA device not detected. Using CPU for torch device.")

jtvae.load_state_dict(torch.load(p.model_loc, map_location=device))

## VAE encoding and decoding on the initial seeding smiles
try:
    design_list = jtvae.smiles_gen(p.seed_smi, p.ndesign)
except KeyError as err:
    print('[KeyError]',err,
          'does not exist in the current JTVAE model vocabulary "%s" \
          (the training set of the model did not contain this structure),'%p.vocab_loc,
          'thus "%s" failed to initialize the model as seeding molecule!'%p.seed_smi)
    exit()

# Check if design generations are ditributed among the top n designs (nseeds)
sub_ndesign = int(p.ndesign)//int(p.nseeds) if p.nseeds > 1 else False

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
    
    ## Report the top design
    top_energy, top_name, top_mol = ranked_poses[0]
    top_smi = top_mol.GetProp('SMILES')
    print("[INFO]: Cycle %s: %s %s kcal/mol"%(j, top_smi, top_energy)+'\t'*6)
    
    ## Generate new design list
    if sub_ndesign:
        design_list = distributed_generator(ranked_poses, p.nseeds, sub_ndesign, jtvae)
    elif p.ensemble > 1:
        top_smi_list = [mol.GetProp('SMILES') for _, _, mol in ranked_poses[:p.ensemble]]
        smi = jtvae.find_ensemble(top_smi_list)
        design_list = jtvae.smiles_gen(smi, p.ndesign)
    else:
        design_list = single_generator(ranked_poses, p.ndesign, jtvae)
    
    
print("\n", p.ncycle, "cycles of design finished. Starting post-processing.")
# Create post-process working directory
postproc_wd = os.path.join(wd, "All_Designs_Processed")
os.makedirs(postproc_wd)
# Extract all ranked designs from each cycle and combine in one sdf file
allmols, bestmols = combine_designs(wd, postproc_wd)
# Create pandas dataframe for summary
allscores, _ = mkdf(allmols, bestmols, postproc_wd)
# Make LSH Forest 
lf = LSH_Convert(allmols, postproc_wd, num_workers = os.cpu_count()-1)
# Get LSH Tree Coords
x, y, s, t = tree_coords(lf, 
                         node_size = float(eval(p.node_size)), 
                         k = int(p.k), 
                         mmm_rps = int(p.mmm_repeats))

# Save coords
with open(os.path.join(postproc_wd,"coords.pickle"),'wb') as f:
    pickle.dump((x,y,s,t),f)
# Create tmap on faerun
f = df_to_faerun(allscores,x,y,s,t)
f.plot("SampleDock"+'_space', path = postproc_wd, # name and path of the .html file
       template="smiles")
with open(os.path.join(postproc_wd,'SampleDock.faerun'), 'wb') as handle:
    pickle.dump(f.create_python_data(), handle, protocol=pickle.HIGHEST_PROTOCOL)