import sys
from sampledock import Vocab, JTNNVAE, hyperparam_loader
import torch

class generator(JTNNVAE):
    def __init__(self, params):
        self.params = params
        vocab = Vocab(self.params.vocab)
        super.__init__(vocab, 
                       self.params.hidden_size, 
                       self.params.latent_size, 
                       self.params.depthT, 
                       self.params.depthG)
        
        if torch.cuda.is_available():
            device = torch.device("cuda")
            print('Using CUDA device:', torch.cuda.get_device_name(torch.cuda.current_device()))
        else:
            device = torch.device("cpu")
        print("CUDA device not detected. Using CPU for torch device.")
        self.load_state_dict(torch.load(self.params.model_loc, map_location=device))
        self.params.sub_ndesign = int(self.params.ndesign)//int(self.params.nseeds) \
        if self.params.nseeds > 1 else False
        self.params.wd = create_wd(self.params.output, self.params.receptor_name)
        
    def generate(self, ranked_poses, ndesign):
        '''
        Generate next cycle of designs based on the sorted poses
        arg:
            - ranked_poses:tuple (energy:float, design:str, best_pose:rdkit.Mol)
            - ndesign:int Number of designs to be generated

        return design_list:list(SMILES:str)
        '''
        for energy, name, mol in ranked_poses:
            smi = mol.GetProp('SMILES')
            try:
                print('[INFO]: Generating new designs \t', end = '\r')
                sys.stdout.flush()
                # get new design list for the nex cycle
                design_list = jtvae.smiles_gen(smi, ndesign)
            # This is due to difference in parsing of SMILES (especially rings)
            ## TODO: Convert sampledock to OOP structure and use the vectors directly
            except KeyError as key:
                print('[KeyError]',key,'is not part of the vocabulary (the model was not trained with this scaffold)')
                continue
            # if there are offspring designs, break the loop
            if len(design_list) != 0: 
                break 
            # go to the next candidate if the current one does not give any return
            else: 
                print('Current design (%s) has no offspring; trying the next one \r'%name)

        return design_list # return a list of SMILES
    
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