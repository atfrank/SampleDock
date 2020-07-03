# Prepare pocket_docking.prm file for rDock
import subprocess
import sys
import os

def mk_dir(parent='.',target='NewTarget'):
    path = os.path.abspath(os.path.join(parent,target))
    i = 0
    while os.path.exists(path):
        i-=1
        path = os.path.abspath(os.path.join(parent,target))+str(i)
    os.makedirs(path)
    print('Directory made: '+path)
    return path

def fill_param(receptor,ligand,recep_name):
    param = ['RBT_PARAMETER_FILE_V1.00\n',
             'TITLE  DOCKING AGAINST %s RECEPTOR\n'%recep_name,
             'RECEPTOR_FILE %s\n'%receptor,
             '\n',
             '##################################################################\n',
             '### CAVITY DEFINITION: REFERENCE LIGAND METHOD\n',
             '##################################################################\n',
             'SECTION MAPPER\n',
             '    SITE_MAPPER RbtLigandSiteMapper\n',
             '    REF_MOL %s\n'%ligand,
             '    RADIUS 7.0\n',
             '    SMALL_SPHERE 1.5\n',
             '    MIN_VOLUME 100\n',
             '    MAX_CAVITIES 1\n',
             '    VOL_INCR 0.0\n',
             '    GRIDSTEP 0.5\n',
             'END_SECTION\n',
             '\n',
             '#################################\n',
             '#CAVITY RESTRAINT PENALTY\n',
             '#################################\n',
             'SECTION CAVITY\n',
             '    SCORING_FUNCTION RbtCavityGridSF\n',
             '    WEIGHT 1.0\n',
             'END_SECTION']
    return param

def wrt_prm(parameter,filename='pocket_docking.prm'):
    with open(filename,'w') as f:
        f.writelines(line for line in parameter)
    f.close()
    print('\nDocking Parameters File Saved at: \n'+filename+'\n')
    sys.stdout.flush()
    
def prep_prm(receptor,ligand,recep_name,target_dir):
    
    receptor = os.path.abspath(receptor)
    if not os.path.exists(receptor): print(receptor+'\nRECEPTOR FILE NOT EXIST!'); return None
    
    ligand = os.path.abspath(ligand)
    if not os.path.exists(ligand): print(ligand+'\nLIGAND FILE NOT EXIST!'); return None
    
    cav_dir = os.path.abspath(target_dir)+'/cavity'
    os.makedirs(cav_dir)
    prmfile = cav_dir+"/pocket_docking.prm"
    param = fill_param(receptor,ligand,recep_name)
    wrt_prm(param,prmfile)
    return prmfile, cav_dir

def create_cav(prmfile):
    # rbcavity must be installed/loaded to execute the cmdline
    cmdline = "rbcavity -was -d -r %s"%prmfile
    proc = subprocess.Popen(cmdline, shell=True)
    proc.wait()
    print('Docking pocket grid created for: \n'+prmfile+'\n')
    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--receptor", help="receptor file path; must be .mol2 format")
    parser.add_argument("-l","--ligand", help="ligand file path; must be .sd format")
    parser.add_argument("-o","--output", help="output directory, default to the current dir",
                       default = os.getcwd())
    parser.add_argument("-t","--target", help="folder name for binding target param files", type=str,
                        default = "new_target")
    a = parser.parse_args()

    target_dir = mk_dir(a.output,a.target)
    prmfile, cav_dir = prep_prm(a.receptor,a.ligand,a.target,target_dir)
    print('\n')
    if prmfile:
        create_cav(prmfile)
    