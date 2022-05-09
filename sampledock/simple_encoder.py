from sampledock import VAEInterface, Design, stockparamloader
import torch
import numpy as np

def encode_smiles_list(smiles_list, vae):
    designs = [Design(smi) for smi in smiles_list]
    for design in designs:
        vae.encode_design(design)
        
    tree_means = torch.cat([design.tree_mean for design in designs]).detach()
    mol_means = torch.cat([design.mol_mean for design in designs]).detach()
    return tree_means, mol_means

def main(smiles_file, model_loc, vocab_loc, outpath):

    model_prms = stockparamloader(model_loc, vocab_loc)
    vae = VAEInterface(model_prms)
    
    with open(smiles_file, 'r') as f:
        smiles_list = [l.strip() for l in f.readlines()]
    print('Number of SMILES to be encoded:', len(smiles_list))
    
    tree_means, mol_means = encode_smiles_list(smiles_list, vae)
    
    np.save(outpath+'_tree_mean', tree_means.numpy())
    np.save(outpath+'_mol_mean', mol_means.numpy())
    print('FINISHED')
    
    
if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--model", help="model location", 
                        default='./models/moses-h450z56/model.iter-400000')
    parser.add_argument("-v", "--vocab", help="vocab location", 
                        default="./models/moses-h450z56/vocab.txt")
    parser.add_argument("-i","--input", help="input smiles file")
    parser.add_argument("-o","--output", help="output .mat files, default to './encoded'",
                    default = './encoded')
    
    a = parser.parse_args()
    main(a.input, a.model, a.vocab, a.output)