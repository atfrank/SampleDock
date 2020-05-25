# Python wrapper for rDock tools
from rdkit import Chem
from rdkit.Chem import AllChem
import os
import subprocess

def dock(ligs, dock_dir, prmfile, docking_prm, npose, prefix = 'docked'):
    # ligs must be a list of file path
    print('Docking in Progress\t', end = '\r')
    procs = []
    for i,lig in enumerate(ligs):
        output = os.path.join(dock_dir,prefix+str(i))
        screenout = os.path.join(dock_dir,prefix+str(i)+'.out')
        cmdline = 'rbdock -i %s -o %s -r %s -p %s -T 1  -n %s > %s'\
                %(lig, output, prmfile, docking_prm, npose, screenout)
        proc = subprocess.Popen(cmdline, shell=True)
        procs.append(proc)

    for proc in procs:
        # makes sure the docking has completed before sorting the score
        proc.wait()
    print('Docking Complete!  \t', end = '\r')

def sort_pose(dock_dir, sort_by, prefix = None):
    # list all pose_org.sd files
    if type(prefix) == str:
        poses_mols = [x for x in os.listdir(dock_dir) 
                      if x.endswith('.sd') and x.startswith(prefix)]
    else:
        poses_mols = [x for x in os.listdir(dock_dir) 
                      if x.endswith('.sd')]
    
    if len(poses_mols) == 0: 
        raise Exception('No .sd file matching the criteria in %s'%dock_dir)
    
    best_poses = []
    for mol_path in poses_mols:
        mol_in = os.path.join(dock_dir,mol_path)
        mol_out = os.path.join(dock_dir,'sorted_'+mol_path)
        cmdline = "sdsort -n -f'%s' %s > %s"%(sort_by, mol_in, mol_out)
        proc = subprocess.Popen(cmdline, shell=True)
        proc.wait()
        # retrieve the best pose mol for each design
        best_pose = Chem.SDMolSupplier(mol_out)[0]
        best_poses.append((float(best_pose.GetProp(sort_by)),best_pose.GetProp('Name'),best_pose))
    print('Docked Poses Sorted       \t', end = '\r')
    return sorted(best_poses)

def save_pose(sorted_poses, save_dir):
    # sort and save the best pose and score for each design
    w = Chem.SDWriter(save_dir+'/ranked_designs.sd')
    f = open(save_dir+'/scores.txt','w')
    
    for score, name, mol in sorted_poses:
        f.write(name+'\t'+str(score)+'\n')
        w.write(mol)
        w.flush()

    w.close()
    f.close()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input", help="directory that contains prepared .sd files")
    parser.add_argument("-o","--output", help="directory for output, default to ./test", 
                        default = "./test")
    parser.add_argument("-r","--prm", help="path to pocket .prm file")
    parser.add_argument("-n","--npose", help="number of poses for docking, default to 100",
                        default = 100)
    parser.add_argument("-s","--sort", help="filter for sdsort, default to SCORE.INTER",
                        default = "SCORE.INTER")
    parser.add_argument("-p","--dprm", help="path to docking prm file, default to dock.prm (no solvation term)",
                        default = "dock.prm")
    parser.add_argument("-x","--prefix", help="rDock output file prefix, default to 'docked'",
                        default = "docked")
    a = parser.parse_args()
    
    a.input = os.path.abspath(a.input)
    a.prm = os.path.abspath(a.prm)
    
    i = 0
    while os.path.exists(a.output):
        i -= 1
        a.output = a.output+str(i)
        
    a.output = os.path.abspath(a.output)
    os.makedirs(a.output)
    
    lig_file_list = [a.input+'/'+x for x in os.listdir(a.input) if x.endswith('.sd')]
    
    dock(lig_file_list, a.output, a.prm, a.dprm, a.npose, a.prefix)
    
    ranked = sort_pose(a.output, a.sort, a.prefix)
    
    save_pose(ranked, a.output)
    
    best_energy, best_mol = ranked[0]
    best_smi = best_mol.GetProp('SMILES')
    print("[INFO]: %s %s kcal/mol"%(best_smi, best_energy))