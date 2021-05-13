from .jtvae import Vocab, JTNNVAE
# Sample and Dock tools
from .SnD import prep_prm
from .SnD import dock, sort_pose, save_pose
from .SnD import hyperparam_loader, create_wd, smiles_to_sdfile
from .SnD import single_generator, distributed_generator
from .SnD import combine_designs, mkdf
from .SnD import LSH_Convert, tree_coords, df_to_faerun

class sampledocker(object):
    def __init__(self, args):
        self.params = hyperparam_loader(args.params)
        self.params.out = args.out
        self.generator = generator(self.params)
        # Create working directory
        self.wd = create_wd(self.params.out, self.params.receptor_name)
        
        ## prepare rdock .prm file and get file path
        self.prmfile, cav_dir = prep_prm(p.receptor_file,p.ligand_file,p.receptor_name,wd)

        ## create pocket
        cmdline = p.cavity_protocol+' %s > %s/create_cavity.out'%(prmfile,cav_dir)
        proc = subprocess.Popen(cmdline, shell=True)
        proc.wait()
        print('Docking pocket grid created')

    def DocknSort(self, smiles_list, npose, sort_by = self.params.sort_by):
        ligs = smiles_to_sdfile(smiles_list, temp_design_dir)
        
        dock(ligand_list, temp_dir, self.prmfile, 
            self.params.docking_prm, npose, self.params.prefix)
        ranked_poses = sort_pose(temp_dir, sort_by, self.params.prefix)
        return ranked_poses
        
    def IterateCycles(self):
        ## VAE encoding and decoding on the initial seeding smiles
        try:
            design_list = self.generator.smiles_gen(self.params.seed_smi, False)
        except KeyError as err:
            raise KeyError(
                err,
                'does not exist in the current JTVAE model vocabulary "%s" \
                (the training set of the model did not contain this structure),'%self.params.vocab_loc,
                'thus "%s" failed to initialize the model as seeding molecule!'%self.params.seed_smi
            )
        ## Main loop: VAE on subsequent returned compounds
        for j in range(self.params.ncycle):

            self.design_dir = os.path.abspath(os.path.join(self.wd,'cycle_%s'%j))
            self.docking_dir = os.path.abspath(os.path.join(self.design_dir, 'docking'))
            os.makedirs(self.docking_dir)

            ## write .sdf file and get ligs file names
            ligs = smiles_to_sdfile(design_list, self.design_dir)

            ## Generate docking scores from .sd files
            dock(ligs, self.docking_dir, self.prmfile, 
                 self.params.docking_prm, self.params.npose, self.params.prefix)
            ranked_poses = sort_pose(self.docking_dir, self.params.sort_by, self.params.prefix) 
            save_pose(ranked_poses, self.design_dir)

            ## Report the top design
            top_energy, top_name, top_mol = ranked_poses[0]
            top_smi = top_mol.GetProp('SMILES')
            print("[INFO]: Cycle %s: %s %s kcal/mol"%(j, top_smi, top_energy)+'\t'*6)

            design_list = self.generator.generate(ranked_poses)
    
    
print("\n", p.ncycle, "cycles of design finished. Starting post-processing.")
# Create post-process working directory
postproc_wd = os.path.join(wd, "All_Designs_Processed")
os.makedirs(postproc_wd)
# Extract all ranked designs from each cycle and combine in one sdf file
allmols, bestmols = combine_designs(wd, postproc_wd)
# Create pandas dataframe for summary
allscores, _ = mkdf(allmols, bestmols, postproc_wd)
# Make LSH Forest 
lf = LSH_Convert(allmols, postproc_wd, num_workers = os.cpu_count()-1)
# Get LSH Tree Coords
x, y, s, t = tree_coords(lf, 
                         node_size = float(eval(p.node_size)), 
                         k = int(p.k), 
                         mmm_rps = int(p.mmm_repeats))

# Save coords
with open(os.path.join(postproc_wd,"coords.pickle"),'wb') as f:
    pickle.dump((x,y,s,t),f)
# Create tmap on faerun
f = df_to_faerun(allscores,x,y,s,t)
f.plot("SampleDock"+'_space', path = postproc_wd, # name and path of the .html file
       template="smiles")
with open(os.path.join(postproc_wd,'SampleDock.faerun'), 'wb') as handle:
    pickle.dump(f.create_python_data(), handle, protocol=pickle.HIGHEST_PROTOCOL)