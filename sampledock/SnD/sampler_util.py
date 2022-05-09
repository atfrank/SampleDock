import torch
import sys
import os

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import rdBase
## Disable rdkit Logs
rdBase.DisableLog('rdApp.error')
from ..jtvae import Vocab, JTNNVAE

import subprocess
from datetime import date

class stockparamloader:
    def __init__(self, model_loc, vocab_loc):
        self.model_loc = model_loc
        vocab_list = [x.strip("\r\n ") for x in open(vocab_loc)]
        self.vocab = Vocab(vocab_list)
        ###### values below are specific to the stock MOSES model #######
        self.hidden_size = 450
        self.latent_size = 56                                   
        self.depthT = 20
        self.depthG = 3

class hyperparamloader:
    def __init__(self,filename=None):
        self.paramfile = os.path.abspath(filename)
        FILE = open(filename)
        print('\n'+'#'*11+' Parameters Loaded as Below '+'#'*11+'\n')
        for i, line in enumerate(FILE):
            if line.strip().startswith("#") or line.isspace(): pass
            else: 
                try:
                    if "#" in line: 
                        line = line.split("#",1)[0]
                    name, value = line.split("=",1)
                    name = name.strip()
                    value = value.strip()
                    if value.isdigit(): value = int(value)
                    setattr(self,name,value)
                    print(name,":",value)
                except:
                    print('\n[ERROR] Failed to load parameter on line %i: \n'%(i+1), line.strip("\r\n "))
                    FILE.close()
                    exit()
        FILE.close()
        print('\n'+'#'*50+'\n')
        self.vocab = [x.strip("\r\n ") for x in open(self.vocab_loc)]
        sys.stdout.flush()
        
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
    sys.stdout.flush()
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
    