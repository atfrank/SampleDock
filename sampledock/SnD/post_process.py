import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import os

from rdkit.Chem import RDConfig # Allow Contrib packages to be used
from rdkit.Chem.Crippen import MolLogP as LogP
from rdkit.Chem.QED import default as QED
from rdkit.Chem.Descriptors import MolWt
import sys
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
from sascorer import calculateScore as SAS

def mkdf(directory,output):
    folders = [x for x in os.listdir(directory) if x.startswith('cycle_')]
    if len(folders) == 0:
        raise Exception('No "cycle_" folder found!')
    scores = pd.DataFrame()
    for i,fd in enumerate(folders):
        df = pd.DataFrame()
        fd_path = os.path.join(directory,fd)
        mols = Chem.SDMolSupplier(fd_path+'/ranked_designs.sd')
        df['Design'] = [m.GetProp('Name') for m in mols]
        df['Cycle'] = i
        df['Score'] = [float(m.GetProp('SCORE.INTER')) for m in mols]
        df['SMILES'] = [m.GetProp('SMILES') for m in mols]
        df['Mol'] = [m for m in mols]
        df['LogP'] = [LogP(m) for m in mols]
        df['QED'] = [QED(m) for m in mols]
        df['MolWt'] = [MolWt(m) for m in mols]
        df['SAS'] = [SAS(m) for m in mols]
        scores = pd.concat([scores,df])

    minscores = scores[scores.index == 0]
    minscores = minscores.sort_values('Score')
    minscores.drop_duplicates('SMILES', inplace = True, keep = 'first')
    scores.to_csv(output+'/all_design.csv')
    minscores.to_csv(output+'/best_designs.csv')
    print("DataFrames Saved!")
    return scores, minscores

def combine_designs(directory, output):
    folders = [x for x in os.listdir(directory) if x.startswith('cycle_')]
    if len(folders) == 0:
        raise Exception('No "cycle_" folder found!')
    mols = []
    best_mols = []
    wa = Chem.SDWriter(output+'/All_Designs.sdf')
    wb = Chem.SDWriter(output+'/Best_Designs.sdf')
    for fd in folders:
        cycle = fd.strip("cycle_")
        sd = directory+'/'+fd+'/ranked_designs.sd'
        if os.path.exists(sd):
            cir_mols = Chem.SDMolSupplier(sd)
            for i, m in enumerate(cir_mols):
                m.SetProp('Cycle',cycle)
                m.SetProp('MolWeight', str(MolWt(m)))
                m.SetProp('LogP', str(LogP(m)))
                m.SetProp('QED', str(QED(m)))
                m.SetProp('SAS', str(SAS(m)))
                mols.append(m)
                wa.write(m)
                if i == 0: 
                    # Select the highest score design in the cycle
                    best_mols.append(m)
                    wb.write(m)
                if int(cycle)%5000 == 0: 
                    wa.flush()
                    wb.flush()
    wa.close()
    wb.close()
    print(len(mols), "total molecules combined from", len(folders),"cycles in\n", directory)
    print(len(best_mols), "selected")
    sys.stdout.flush()
    return mols, best_mols

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="combine and the ranked_designs.sd in each "+
                                     "'cycle_*' folder from Sample and Dock and calculate MolWeight, SAS, LogP, and QED.")
    parser.add_argument("-i","--input", help="input directory that contain folder by cycles")
    parser.add_argument("-o","--outpath", help="output directory for the combined sdf file",
                        default='./')
    a = parser.parse_args()
    directory = os.path.abspath(a.input)
    out = os.path.abspath(a.outpath)
    if not os.path.exists(out): 
        os.makedirs(out)
        print(out, "Made")
    combine_designs(directory, out)
    mkdf(directory, out)
