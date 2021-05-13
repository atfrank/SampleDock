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
        ## Load model weights
        self.load_state_dict(torch.load(self.params.model_loc, map_location=device))
        self.distributed = (self.params.nseeds > 1)
        if self.distributed:
            self.params.ndesign = int(self.params.ndesign)//int(self.params.nseeds)

    ## Gaussian sampling and smiles generation for a given smiles
    def smiles_gen(self, smiles, tree_mean = None, mol_mean = None):
        ## Convert smiles to one-hot encoding (altered function from original code)
        x_tree, x_mol = self.encode_single_smiles(smiles)
        
        ## Encode one-hots to z-mean and log variance.
        if tree_mean == None and mol_mean == None:
            tree_mean = self.T_mean(x_tree)
            mol_mean = self.G_mean(x_mol)
        
        tree_log_var = -torch.abs(self.T_var(x_tree)) 
        mol_log_var = -torch.abs(self.G_var(x_mol))

        smiles_list = []
        vec_list = []
        for i in range(self.ndesign):
            ## generate latent vectors (stochastic)
            z_tree = self.add_noise_z(tree_mean, tree_log_var)
            z_mol = self.add_noise_z(mol_mean, mol_log_var)
            ## decode back to smiles
            smilesout = self.decode(z_tree, z_mol, prob_decode = False)
            ## Check if the smiles already exists
            if smilesout not in smiles_list:
                smiles_list.append(smilesout)
                vec_list.append((z_tree, z_mol))
        return smiles_list, vec_list
        
    def find_ensemble(self, smiles_list):
        z_tree = []
        z_mol = []
        for smi in smiles_list:
            try:
                x_tree, x_mol = self.encode_single_smiles(smi)
            # This is due to difference in parsing of SMILES (especially rings)
            ## TODO: Convert sampledock to OOP structure and use the vectors directly
            except KeyError as key:
                print('[KeyError]',key,'is not part of the vocabulary (the model was not trained with this scaffold)')
                continue
            tree_mean = self.T_mean(x_tree)
            tree_log_var = -torch.abs(self.T_var(x_tree)) 
            mol_mean = self.G_mean(x_mol)
            mol_log_var = -torch.abs(self.G_var(x_mol))
            
            z_tree.append(self.add_noise_z(tree_mean, tree_log_var))
            z_mol.append(self.add_noise_z(mol_mean, mol_log_var))
        
        z_tree = torch.cat(z_tree)
        z_mol = torch.cat(z_mol)
        
        return self.decode(z_tree.mean(0).reshape((1,self.latent_size)),
                            z_mol.mean(0).reshape((1,self.latent_size)),
                            False)
    
    def _generator(self, smi):
        try:
            print('[INFO]: Generating new designs \t', end = '\r')
            sys.stdout.flush()
            # get new design list for the nex cycle
            design_list = self.smiles_gen(smi, self.params.ndesign)
            return design_list # return a list of SMILES

        # This is due to difference in parsing of SMILES (especially rings)
        ## TODO: Convert sampledock to OOP structure and use the vectors directly 
        except KeyError as key:
            print('[KeyError]',key,'is not part of the vocabulary (the model was not trained with this scaffold)')
            
    def generate(self, ranked_poses):
        '''
        Generate next cycle of designs based on the sorted poses
        arg:
            - ranked_poses:tuple (energy:float, design:str, best_pose:rdkit.Mol)

        return design_list:list(SMILES:str)
        '''
        ## Generate new design list
        if self.distributed:
            design_list = []
            for i_seed in range(self.params.nseeds):
                energy, name, mol = ranked_poses[i_seed]
                smi = mol.GetProp('SMILES')
                cur_design_list = self._generator(smi)
                if len(cur_design_list) > 0:
                    design_list.extend(cur_design_list)
                    
        elif self.params.ensemble > 1:
            top_smi_list = [mol.GetProp('SMILES') for _, _, mol in ranked_poses[:self.params.ensemble]]
            smi = self.find_ensemble(top_smi_list)
            design_list = self.smiles_gen(smi)
            
        else:
            for energy, name, mol in ranked_poses:
                smi = mol.GetProp('SMILES')
                design_list = self._generator(smi)
                if len(design_list) > 0: 
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