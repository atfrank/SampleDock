import torch
import sys
import os

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import rdBase
## Disable rdkit Logs
rdBase.DisableLog('rdApp.error')

import subprocess
from datetime import date
# JTVAE tools
from jtvae import Vocab, JTNNVAE
# Sample and Dock tools
from pocket_prepare import prep_prm
from docking import dock, sort_pose, save_pose

class hyperparam_loader(object):
    def __init__(self,filename="hyper.param"):
        self.paramfile = os.path.abspath(filename)
        FILE = open(filename)
        print('\n'+'#'*11+' Parameters Loaded as Below '+'#'*11+'\n')
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
        print('\n'+'#'*50+'\n')
        self.vocab = [x.strip("\r\n ") for x in open(self.vocab_loc)]
        
def create_wd(parent_dir,target_name):
    
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
    print("\nNew Directory Made:"+directory)
    return directory

def smiles_to_sdfile(smiles_list, dsgn_dir):
    lig_file_names = []
    for i, smi in enumerate(smiles_list):
        name = 'design_'+str(i)
        output = dsgn_dir+'/'+name+'.sd'
        m2 = Chem.MolFromSmiles(smi)
        AllChem.Compute2DCoords(m2)
        m2.SetProp("_Name", name)
        m3 = Chem.AddHs(m2)
        AllChem.EmbedMolecule(m3,AllChem.ETKDG())
        m3.SetProp("SMILES", smi)
        w = Chem.SDWriter(output)
        w.write(m3)
        w.flush()
        lig_file_names.append(output)
    w.close()
    return lig_file_names
        
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-p","--params", help="parameter file for SnD", 
                        default = "./hyper.param")
    parser.add_argument("-o","--output", help="directory for output, default to ../SnD_designs", 
                        default = "../SnD_designs")
    a = parser.parse_args()
    
    # Load hyper parameters
    p = hyperparam_loader(a.params)
    
    ## Load Stock JTNN VAE Model
    vocab = Vocab(p.vocab)
    jtvae = JTNNVAE(vocab, p.hidden_size, p.latent_size, p.depthT, p.depthG)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    jtvae.load_state_dict(torch.load(p.model_loc, map_location=device))
    
    ## VAE encoding and decoding on the initial seeding smiles
    try:
        design_list = jtvae.smiles_gen(p.seed_smi, p.ndesign)
    except KeyError as err:
        print('KeyError:',err,
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
        try: os.makedirs(docking_dir)
        except FileExistsError: print(docking_dir,'Overwritten')

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
                design_list = jtvae.smiles_gen(smi, p.ndesign)
            # go to the second best candidate if the best does not give any return
            except KeyError as err:
                print('KeyError:',err,'is not part of the vocabulary')
                continue
            
            if len(design_list) != 0: 
                break 
                
            else: 
                print('Current design (%s) has no offspring; trying the next one \r'%name)
                
        print("[INFO]: Cycle %s: %s %s kcal/mol"%(j, smi, energy)+'\t'*6)
