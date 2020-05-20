import torch
import torch.nn as nn
import sys
sys.path.append('/home/ziqiaoxu/Sample_Dock_DnD/jtvae/')
import os
import argparse

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import rdBase
## Disable rdkit Logs
rdBase.DisableLog('rdApp.error')

import subprocess
from datetime import date

from torch.autograd import Variable
from nnutils import create_var
from mol_tree import Vocab
from jtnn_vae import JTNNVAE
from pocket_prepare import prep_prm

class param_loader(object):
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
    os.makedirs(directory+'/design/')
    print("New Directory Made:"+directory)
    return directory

def z_vecs(x_mean, x_log_var):
    
    epsilon = create_var(torch.randn_like(x_mean))
    z_vecs = x_mean + torch.exp(x_log_var / 2) * epsilon
    return z_vecs

def smiles_gen(smiles,ndesigns):
    ## Convert smiles to one-hot encoding (altered function from original code)
    x_tree, x_mol = jtnn.encode_from_smiles_xs(smiles)
    ## Encode one-hots to z-mean and log var. Following Mueller et al.
    tree_mean = jtnn.T_mean(x_tree)
    tree_log_var = -torch.abs(jtnn.T_var(x_tree)) 
    mol_mean = jtnn.G_mean(x_mol)
    mol_log_var = -torch.abs(jtnn.G_var(x_mol))

    smiles_list = []
    for i in range(ndesigns):
        ## generate latent vectors (stochastic)
        z_tree = z_vecs(tree_mean, tree_log_var)
        z_mol = z_vecs(mol_mean, mol_log_var)
        ## decode back to smiles
        smilesout = jtnn.decode(z_tree,z_mol,False)
        ## Check if the smiles already exists
        if smilesout not in smiles_list:
            smiles_list.append(smilesout)    
    return smiles_list

def smiles_to_sdfile(smiles_list,directory):
    for i, x in enumerate(smiles_list):
        name = 'design_'+str(i)
        output = '%s/design/'%directory+name+'.sd'
        m2 = Chem.MolFromSmiles(x)
        AllChem.Compute2DCoords(m2)
        m2.SetProp("_Name", name)
        m3 = Chem.AddHs(m2)
        AllChem.EmbedMolecule(m3,AllChem.ETKDG())
        w = Chem.SDWriter(output)
        w.write(m3)
        w.flush()
        
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-p","--params", help="parameter file for SnD", 
                        default = "./hyper.param")
    parser.add_argument("-o","--output", help="directory for output,default to cir", 
                        default = "../designs")
    a = parser.parse_args()
    
    # Load hyper parameters
    p = param_loader(a.params)
    
    # Create working directory
    wd = create_dirs(a.output,p.receptor_name)
    
    ## Load Stock JTNN VAE Model
    vocab = Vocab(p.vocab)
    jtnn = JTNNVAE(vocab, p.hidden_size, p.latent_size, p.depthT, p.depthG)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    jtnn.load_state_dict(torch.load(p.model_loc, map_location=device))
    
    ## prepare rdock .prm file
    prmfile = prep_prm(p.receptor_file,p.ligand_file,p.receptor_name,wd)
    
    ## create pocket
    cmdline = p.cavity_protocol+" "+prmfile
    proc = subprocess.Popen(cmdline, shell=True)
    proc.wait()
    print('Docking pocket grid created for: \n'+prmfile+'\n')
    
    ## Parameters for loops
    cycle = p.ncycle  # design generation cycles
    ndesigns = p.ndesign  # maximum number of designs should be generated

    score_file = "%s/design/scores.txt"%wd

    ## VAE encoding and decoding on the first input smiles
    design_list = smiles_gen(p.seed_smi, ndesigns)
    
    ## Main loop: VAE on subsequent returned compounds
    for j in range(cycle):

        ## write .sdf file
        smiles_to_sdfile(design_list,wd)

        ## Generate docking scores from .sdf files
        cmdline = "/usr/bin/bash vae_docking.sh %s %s/design %s"%(len(design_list),wd,prmfile)
        proc = subprocess.Popen(cmdline, shell=True)
        proc.wait()

        # process scores and set best_design
        scores = pd.read_csv(score_file, sep = " ", names = ['design', 'score'])
        scores.sort_values(by = 'score', inplace = True)

        ## Select the best score 
        best_energy = scores.score.values[0]
        best_design_idx = scores.design.values[0]
        best_design = design_list[best_design_idx]

        ## Encode and decode by VAE and update design_list
        design_list = smiles_gen(best_design, ndesigns)

        ## Copy working files
        cmdline = "cp -r %s/design/ %s/design_%s/"%(wd,wd,j)
        proc = subprocess.Popen(cmdline, shell=True)
        proc.wait()

        ## Clean working dir content
        cmdline = "rm -r %s/design/*"%wd
        proc = subprocess.Popen(cmdline, shell=True)
        proc.wait()

        print("[INFO]: Cycle %s: %s %s kcal/mol"%(j, best_design, best_energy))


    cmdline = "bash combine_sd_in_sampler.sh %s"%wd
    proc = subprocess.Popen(cmdline, shell=True)
    proc.wait()