import sys
from sampledock import JTNNVAE
from sampledock.jtvae.mol_tree import MolTree
from sampledock.jtvae.nnutils import create_var
from sampledock.jtvae.datautils import tensorize
import torch
import numpy as np
from random import random

class Design:
    """
    Data container for designs
    Holding smiles, latent vectors, log_var, and score
    """
    def __init__(self, smiles):
        self.smiles = smiles
        self.tree_mean = None
        self.mol_mean = None
        self.tree_log_var = None 
        self.mol_log_var = None
        self.score = None

class VAEInterface(JTNNVAE):
    def __init__(self, model_prms):
        super().__init__(model_prms.vocab, model_prms.hidden_size, 
                         model_prms.latent_size, model_prms.depthT, model_prms.depthG)
        device = torch.device('cuda' if torch.cuda.is_available() else "cpu")
        super().load_state_dict(torch.load(model_prms.model_loc, map_location=device))
    
    def encode_design(self, design: Design) -> None:
        tree_batch = [MolTree(design.smiles)]
        _, jtenc_holder, mpn_holder = tensorize(tree_batch, self.vocab, assm=False)
        tree_x_vecs, _, mol_x_vecs = self.encode(jtenc_holder, mpn_holder)
        design.tree_mean = self.T_mean(tree_x_vecs)
        design.mol_mean = self.G_mean(mol_x_vecs)
        design.tree_log_var = -torch.abs(self.T_var(tree_x_vecs))  
        design.mol_log_var = -torch.abs(self.G_var(mol_x_vecs)) 

    def perturb_decode(self, design: Design, sigma: int = 2) -> Design:
        tree_epsilon = create_var(torch.randn_like(design.tree_mean))
        mol_epsilon = create_var(torch.randn_like(design.mol_mean))
        # tree_mean = design.tree_mean + torch.exp(design.tree_log_var / 2) * tree_epsilon * sigma
        # mol_mean = design.mol_mean + torch.exp(design.mol_log_var / 2) * mol_epsilon
        tree_mean = design.tree_mean +  tree_epsilon * sigma
        mol_mean = design.mol_mean + mol_epsilon
        new_design =Design(self.decode(tree_mean, mol_mean,prob_decode = False))
        new_design.tree_mean = tree_mean
        new_design.mol_mean = mol_mean
        return new_design

class Lineage:
    """
    Lineage operates on a single design 
    and generates nsamples of neigbor molecules each iteration
    It is initialized with a seeding design, samples the neighborhood iteratively,
    and updates the seeding design.
    The best design in the neighborhood determined by oracle is compared to the seed,
    and the better one becomes the seeding molecule.
    acc_rate specify the chance of accepting low scoring neighbor design as seed.
    """
    def __init__(self, nsamples, vae: VAEInterface, oracle, acc_rate,  sigma, sort_descending):
        self.nsamples = nsamples 
        self.vae = vae
        self.oracle = oracle
        self.acc_rate = acc_rate
        self.design = self.initial_sample()
        self.sort_descending = sort_descending
        self.sigma = sigma

    def assign_score(self, design: Design):
        design.score = self.oracle(design.smiles)

    def initial_sample(self) -> Design:
        # Random Sampling in the latent space
        design = Design(self.vae.sample_prior())
        try:
            # Encoding may fail due to key error
            self.vae.encode_design(design)
            self.assign_score(design) 
            return design

        except KeyError: 
            # resample if encoding failed
            # This is due to difference in parsing of SMILES (especially rings)
            return self.initial_sample()

    def perturb_generate(self) -> Design:
        '''
        Perturb the latent vectors and decode to new smiles
        Update the smiles and regenerate latent vectors based on the new smiles
        '''
        new_design = self.vae.perturb_decode(self.design, self.sigma)
        try:
            # self.vae.encode_design(new_design)
            self.assign_score(new_design)
            return new_design
        
        except KeyError as key:
            print(", ".join(key.args)+' is not part of the vocabulary')
            return self.perturb_generate()

    def produce_samples(self):
        self.samples = [self.perturb_generate() for i in range(self.nsamples)]
        sorted(self.samples, key = lambda design: design.score[0], reverse=self.sort_descending)
        best_sample = self.samples[0]
        if best_sample.score[0] >= self.design.score[0] or random() < self.acc_rate: 
            self.design = best_sample
    
class Generator:
    def __init__(self, batch_size, nsamples, vae : VAEInterface,
                 oracle, acc_rate, steps, sigma, sort_descending=True,
                 ):
        
        self.cur_step = 0
        self.max_steps = steps
        self.lineages = [
            Lineage(nsamples, vae, oracle, 
                    acc_rate, sigma, sort_descending)\
            for i in range(batch_size)
        ]

    def update_batch(self):
        for lineage in self.lineages:
            lineage.produce_samples()

    def __iter__(self):
        return self

    def __next__(self):
        if self.cur_step < self.max_steps:
            self.update_batch()
            self.cur_step += 1
            return self.lineages
        else:
            raise StopIteration

if __name__ == "__main__":
    from sampledock.SnD import stockparamloader
    from rdkit.Chem.Descriptors import MolWt
    from rdkit import Chem
    import argparse
    import time
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--model", help="model location", 
                        default='./models/moses-h450z56/model.iter-400000')
    parser.add_argument("-v", "--vocab", help="vocab location", 
                        default="./models/moses-h450z56/vocab.txt")
    parser.add_argument("-i","--input", help="input smiles",
                        default='c1ccccc1')
    parser.add_argument("-n","--ndesign", help="number of designs to be generated, default to 10",
                    default = 10, type=int)
    parser.add_argument("-o","--output", help="output .smi file, default to './generated.smi'",
                    default = './generated.smi')
    parser.add_argument("-s", "--sigma", type = int, 
                        help = "sampling width in the latent space, default = 2",
                        defualt = 2)
    a = parser.parse_args()
    torch.manual_seed(42)
    model_prms = stockparamloader(a.model,a.vocab)

    def molwt(smiles):
        return MolWt(Chem.MolFromSmiles(smiles))
    
    vae = VAEInterface(model_prms)

    generator = Generator(vae = vae ,batch_size = 10, 
                        nsamples=10, steps = 10, oracle=molwt,
                        sort_descending=False, acc_rate=0.05)
    start_time = time.time()
    for lineages in generator:
        print('Step:', generator.cur_step)
        for lineage in lineages:
            print(round(lineage.design.score ,2),'\t',lineage.design.smiles)
        print('Time Used:', time.time()-start_time)
        print('***'*10)
        start_time = time.time()