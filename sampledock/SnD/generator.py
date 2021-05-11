import sys
from sampledock import Vocab, JTNNVAE
import torch

def jtvae_loader(params):
    vocab = Vocab(params.vocab)
    jtvae = JTNNVAE(vocab, params.hidden_size, params.latent_size, params.depthT, params.depthG)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    jtvae.load_state_dict(torch.load(params.model_loc, map_location=device))
    return jtvae

def single_generator(ranked_poses, ndesign, jtvae):
    '''
    Generate next cycle of designs based on the top scoring pose
    arg:
        - ranked_poses:tuple (energy:float, design:str, best_pose:rdkit.Mol)
        - ndesign:int Number of designs to be generated
    
    return design_list:list(SMILES:str)
    '''
    for energy, name, mol in ranked_poses:
        smi = mol.GetProp('SMILES')
        design_list = _generator(smi, ndesign, jtvae)
        if len(design_list) > 0: 
            break 
        # go to the next candidate if the current one does not give any return
        else: 
            print('Current design (%s) has no offspring; trying the next one \r'%name)

    return design_list # return a list of SMILES

def _generator(smi, ndesign, jtvae):
    try:
        print('[INFO]: Generating new designs \t', end = '\r')
        sys.stdout.flush()
        # get new design list for the nex cycle
        design_list = jtvae.smiles_gen(smi, ndesign)
        return design_list # return a list of SMILES
    
    # This is due to difference in parsing of SMILES (especially rings)
    ## TODO: Convert sampledock to OOP structure and use the vectors directly 
    except KeyError as key:
        print('[KeyError]',key,'is not part of the vocabulary (the model was not trained with this scaffold)')

    

def distributed_generator(ranked_poses, nseeds, sub_ndesign, jtvae):
    '''
    Generate next cycle of designs based on the top n scoring poses
    arg:
        - ranked_poses:tuple (energy:float, design:str, best_pose:rdkit.Mol)
        - nseeds:int Number of the top designs being used
        - sub_ndesign:int Number of designs to be generated for each seed
    
    return design_list:list(SMILES:str)
    '''
    design_list = []
    for i_seed in range(nseeds):
        energy, name, mol = ranked_poses[i_seed]
        smi = mol.GetProp('SMILES')
        cur_design_list = _generator(smi, sub_ndesign, jtvae)
        if len(cur_design_list) > 0:
            design_list.extend(cur_design_list)
            
    return design_list
        
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