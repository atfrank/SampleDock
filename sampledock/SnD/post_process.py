# Post processing script for sample and dock generated molecules

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import os
from multiprocessing import Pool
from itertools import repeat

from rdkit.Chem.PropertyMol import PropertyMol # Allow pickle on mol props for multiprocessing
from rdkit.Chem import RDConfig # Allow Contrib packages to be used
from rdkit.Chem.Crippen import MolLogP as LogP # Lipophilicity
from rdkit.Chem.QED import default as QED # Quantitiative Estimate of Drug-likeness
from rdkit.Chem.Descriptors import MolWt # Mol Weight
import sys
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
# add path for rdkit Contrib packages
from sascorer import calculateScore as SAS # Sythetic Accessiblilty Score

# Function for calculate mol properties for sd files in each folder for multiprocessing
def process_by_folder(fd, inpath):
    cycle = fd.strip("cycle_")
    sd = inpath+'/'+fd+'/ranked_designs.sd'
    if os.path.exists(sd):
        cir_mols = [PropertyMol(m) for m in Chem.SDMolSupplier(sd)]
        for i,m in enumerate(cir_mols):
            # Calculate properties for each mol
            m.SetProp('Cycle',cycle)
            m.SetProp('MolWeight', str(MolWt(m)))
            m.SetProp('LogP', str(LogP(m)))
            m.SetProp('QED', str(QED(m)))
            m.SetProp('SAS', str(SAS(m)))
            if i == 0: 
                # Select the highest score design in the cycle
                best_mol = m
    return cir_mols, best_mol

# calculated mol properties from each cycle and combine mols in one sdf file
def combine_designs(inpath, outpath):
    # list the folders in the directory for all cycles
    folders = [x for x in os.listdir(inpath) if x.startswith('cycle_')]
    # sort folder name
    folders.sort(key=lambda x: int(x.strip('cycle_')))

    if len(folders) == 0:
        raise Exception('No "cycle_" folder found!')
    
    # Multiprocessing
    with Pool(processes = os.cpu_count()-1) as pool:
        results = pool.starmap(process_by_folder, zip(folders, repeat(inpath)))

    # Retrieve results
    mol_lists, best_mols = zip(*results)
    # Create the list of all mols
    all_mols = []
    for l in mol_lists:
        all_mols.extend(l)
    # Convert tuple to list
    best_mols = list(best_mols)

    print(len(all_mols), "total molecules combined from", len(folders),"cycles in\n", inpath)
    print(len(best_mols), "best designs extracted.\n")
    sys.stdout.flush()

    # Save as sdf
    with open(outpath+'/All_Designs.sdf','w') as outfile:
        w = Chem.SDWriter(outfile)
        for m in all_mols:
            w.write(m)
        w.close()

    with open(outpath+'/Best_Designs.sdf','w') as outfile:
        w = Chem.SDWriter(outfile)
        for m in best_mols:
            w.write(m)
        w.close()
    print('Mols saved!')
    sys.stdout.flush()

    return all_mols, best_mols

# Create dataframe with all the properties
def create_df(mol_list):
    df = pd.DataFrame()

    df['Design'] = [m.GetProp('Name') for m in mol_list]
    df['Cycle'] = [int(m.GetProp('Cycle')) for m in mol_list]
    df['Score'] = [float(m.GetProp('SCORE.INTER')) for m in mol_list]
    df['SMILES'] = [m.GetProp('SMILES') for m in mol_list]
    df['Mol'] = [m for m in mol_list]
    df['LogP'] = [float(m.GetProp('LogP')) for m in mol_list]
    df['QED'] = [float(m.GetProp('QED')) for m in mol_list]
    df['MolWt'] = [float(m.GetProp('MolWeight')) for m in mol_list]
    df['SAS'] = [float(m.GetProp('SAS')) for m in mol_list]

    return df

def mkdf(all_mols, best_mols, outpath):
    # Create dataframe from the lists
    allscores = create_df(all_mols)
    minscores = create_df(best_mols)

    # sort the dataframe based on docking scores
    sortedscores = minscores.sort_values('Score')
    # Drop dulicated entries
    sortedscores.drop_duplicates('SMILES', inplace = True, keep = 'first')

    # Save as csv
    allscores.drop(columns=['Mol']).to_csv(outpath+'/allscores.csv', index = False)
    sortedscores.drop(columns=['Mol']).to_csv(outpath+'/sortedscores.csv', index = False)
    print('Dataframes saved!')
    sys.stdout.flush()
    return allscores, minscores

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="combine and the ranked_designs.sd in each "+
                                     "'cycle_*' folder from Sample and Dock and calculate MolWeight, SAS, LogP, and QED.")
    parser.add_argument("-i","--input", help="input directory that contain folder by cycles")
    parser.add_argument("-o","--outpath", help="output directory for the combined sdf file,"+\
                        "default to ./processed_data")
    a = parser.parse_args()
    inpath = os.path.abspath(a.input)

    if a.outpath:
        outpath = os.path.abspath(a.outpath)
    else: outpath = inpath+"/All_Designs_Processed/"
    
    if not os.path.exists(outpath): 
        os.makedirs(outpath)
        print("Directory Made:")
        print(outpath)
        sys.stdout.flush()
    allmols, bestmols = combine_designs(inpath, outpath)
    mkdf(allmols, bestmols, outpath)
