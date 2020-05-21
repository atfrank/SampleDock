import torch
import torch.nn as nn
import sys
sys.path.append('/home/ziqiaoxu/Sample_Dock_DnD/jtvae/')
import os
import argparse

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import rdBase
## Disable rdkit Logs
rdBase.DisableLog('rdApp.error')

import subprocess
from datetime import date
# JTVAE tools
from torch.autograd import Variable
from nnutils import create_var
from mol_tree import Vocab
from jtnn_vae import JTNNVAE
# Sample and Dock tools
from pocket_prepare import prep_prm
from docking import smiles_to_sdfile, dock, sort_pose, save_pose

class hyperparam_loader(object):
    def __init__(self,filename="SnD.param"):
        FILE = open(filename)
        print('Parameters are loaded as below\n')
        for l in FILE:
            if l.startswith("#") or l.isspace(): pass
            else: 
                if "#" in l: l = l.split("#",1)[0]
                name, value = l.split("=",1)
                name = name.strip("\r\n ")
                value = value.strip("\r\n ")
                if value.isdigit(): value = int(value)
                setattr(self,name,value)
                print(name,":",value)
        FILE.close()
        self.vocab = [x.strip("\r\n ") for x in open(self.vocab_loc)]
        
def create_dirs(parent_dir,target_name):
    ## Create working directory marked by date
    td = date.today().strftime("%b%d")
    directory = os.path.join(parent_dir,"SnD-%s-%s"%(target_name,td))
    ## Renaming to avoid overwriting existing data
    i = 0
    while os.path.exists(directory):
        i -= 1
        directory = os.path.join(parent_dir,"SnD-%s-%s"%(target_name,td)+str(i))
    directory = os.path.abspath(directory)
    os.makedirs(directory)
    print("New Directory Made:"+directory)
    return directory


        
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-p","--params", help="parameter file for SnD", 
                        default = "./hyper.param")
    parser.add_argument("-o","--output", help="directory for output, default to ../test", 
                        default = "../test")
    a = parser.parse_args()
    
    # Load hyper parameters
    p = hyperparam_loader(a.params)
    
    # Create working directory
    wd = create_dirs(a.output,p.receptor_name)
    
    ## Load Stock JTNN VAE Model
    vocab = Vocab(p.vocab)
    jtvae = JTNNVAE(vocab, p.hidden_size, p.latent_size, p.depthT, p.depthG)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    jtvae.load_state_dict(torch.load(p.model_loc, map_location=device))
    
    ## prepare rdock .prm file and get file path
    prmfile = prep_prm(p.receptor_file,p.ligand_file,p.receptor_name,wd)
    
    ## create pocket
    cmdline = p.cavity_protocol+" "+prmfile
    proc = subprocess.Popen(cmdline, shell=True)
    proc.wait()
    print('Docking pocket grid created for: \n'+prmfile+'\n')

    ## VAE encoding and decoding on the first input smiles
    design_list = jtvae.smiles_gen(p.seed_smi, p.ndesign)

    ## Main loop: VAE on subsequent returned compounds
    for j in range(p.ncycle):

        design_dir = os.path.abspath(os.path.join(wd,'cycle_%s'%j))
        docking_dir = os.path.abspath(os.path.join(design_dir, 'docking'))
        try: os.makedirs(docking_dir)
        except FileExistsError: print(docking_dir,'Overwritten')

        ## write .sdf file and get ligs file names
        ligs = smiles_to_sdfile(design_list,design_dir)

        ## Generate docking scores from .sd files
        dock(ligs, docking_dir, prmfile, p.docking_prm, p.npose)
        ranked_poses = sort_pose(docking_dir, p.sort_by)
        save_pose(ranked_poses, design_dir)

        ## determine the winner
        best_energy, best_mol = ranked_poses[0]
        best_smi = best_mol.GetProp('SMILES')
        print("[INFO]: Cycle %s: %s %s kcal/mol"%(j, best_smi, best_energy))

        design_list = jtvae.smiles_gen(best_smi, p.ndesign)
