import sys
import os
sys.path.append('/home/ziqiaoxu/Sample_Dock_DnD/jtvae/')

from mol_tree import Vocab
from jtnn_vae import JTNNVAE
import torch

def jtvae_loader(params):
    vocab = Vocab(params.vocab)
    jtvae = JTNNVAE(vocab, params.hidden_size, params.latent_size, params.depthT, params.depthG)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    jtvae.load_state_dict(torch.load(params.model_loc, map_location=device))
    return jtvae

if __name__ == "__main__":
    from sampler import hyperparam_loader
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--params", help="hyperparameter file path")
    parser.add_argument("-i","--input", help="input smiles")
    parser.add_argument("-n","--ndesign", help="number of designs to be generated, default to 10",
                       default = 10)
    parser.add_argument("-o","--output", help="output .smi file, default to './generated.smi'",
                       default = './generated.smi')
    a = parser.parse_args()
    
    jtvae = jtvae_loader(hyperparam_loader(a.params))
    design_list = jtvae.smiles_gen(a.input, a.ndesign)
    print('Seeding SMILES:',a.input)
    with open(a.output,'w') as f:
        for smi in design_list:
            print(smi)
            f.write(smi+'\n')
    f.close()