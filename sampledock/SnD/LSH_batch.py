# Batch Processing SMILES strings for LSH Forest data
# .smi file with one smiles string per line is required
# This is part of SampleDock package
#
# COMMERCIAL LICENSE: https://github.com/atfrank/SampleDock/blob/master/LICENSE_COMMERICAL
# NON-COMMERCIAL LICENSE: https://github.com/atfrank/SampleDock/blob/master/LICENSE_NON-COMMERICAL
# Frank Lab 2020-2021

import sys
import os
import pickle
#from multiprocessing import Pool
from tqdm.contrib.concurrent import process_map
from collections import namedtuple

import tmap as tm
from mhfp.encoder import MHFPEncoder
from rdkit.Chem import RDConfig
from rdkit.Chem import SDMolSupplier, MolFromSmiles, MolToSmiles
from rdkit.Chem.Crippen import MolLogP as LogP
from rdkit.Chem.QED import default as QED
from rdkit.Chem.Descriptors import MolWt
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
from sascorer import calculateScore as SAS

def smi_to_mol(smi):
    # Wrapper function for proper multiprocessing on converting SMILES to mol
    # Chem.MolFromSmiles is Boost.Python.function and cannot be pickled for multiprocessing
    mol = MolFromSmiles(smi)
    return mol

def file_to_mols(filepath):
    if filepath.endswith('.smi'):
        print('Converting SMILES to list of Mols')
        sys.stdout.flush()
        with open(filepath) as infile:
            smiles_list = [line.rstrip() for line in infile.readlines()]
        # Multiprocessing with all available threads
        #with Pool(processes = os.cpu_count()) as pool:
            #mols = pool.map(smi_to_mol, smiles_list)
            
        mols = process_map(smi_to_mol, smiles_list, chunksize = 100, 
                          max_workers = a.worker)
        
        mols = [m for m in mols if m]
                
    elif filepath.endswith('.sd') or filepath.endswith('.sdf'):
        mols = [mol for mol in SDMolSupplier(filepath) if mol]
        
    else: 
        raise Exception('Invalid file: {}\n'.format(filepath)+
             '.smi, .sd, or .sdf extension is expected')
        
    return mols

def cal_fp_props(mol):
    # Wrapper function for multiprocessing
    mol_props = Props(MolToSmiles(mol),MolWt(mol),LogP(mol),QED(mol),SAS(mol))
    fp = enc.encode_mol(mol)
    return fp, mol_props

def single_convert(mol_list, outpath, props_named_tuple = None):
    # Multiprocessing with all available threads
    #with Pool(processes = os.cpu_count()) as pool:
        #fps_and_props = pool.map(cal_fp_props, mol_list)
        
    fps_and_props = process_map(cal_fp_props, mol_list, 
                                chunksize = 100, 
                                max_workers = a.worker)
    
    # Unpack the returned results
    fps, props = zip(*fps_and_props)

    # Turn props into list of named tuples
    if props_named_tuple:
        props = [props_named_tuple(*p) for p in props]

    with open(os.path.join(outpath,"props.pickle"), "wb") as pf:
        pickle.dump(props,pf)
    with open(os.path.join(outpath,"fps.pickle"), "wb+") as f:
        pickle.dump(fps,f)
    
    return fps, props

def batch_divider(length,batch_size):
    # For better reporting progress and including the remainders
    # Define the named tuple for indice
    Batch = namedtuple('Batch',['batch_num','start','stop'])

    # Create an inclusive list of the batch start and stop indice marks
    idx = [*range(0,length,batch_size),length]

    # List of named tuples
    batches = [Batch(*b) for b in zip(range(1,len(idx)),idx[:-1],idx[1:])]
    return batches

def batch_convert(mol_list, batch_size, outpath, props_named_tuple = None):
    # Divide batches
    batches = batch_divider(len(mol_list), batch_size)
    
    propfile = open(os.path.join(outpath,"props.pickle"), "wb")
    fpfile = open(os.path.join(outpath,"fps.pickle"), "wb")
    
    for batch_num, start, stop in batches:
        # Reporting progress
        print('Current batch: {} of {} \t'.format(batch_num,len(batches)),
             end = '\r')
        sys.stdout.flush()
        batch_mol_list = mol_list[start:stop]
        # Multiprocessing with all available threads
        #with Pool(processes = os.cpu_count()) as pool:
            #fps_and_props = pool.map(cal_fp_props, batch_mol_list)
        
        # Use process map for a progress bar
        fps_and_props = process_map(cal_fp_props, batch_mol_list, 
                                    chunksize = 100, 
                                    max_workers = a.worker)
        
        # Unpack the returned results
        fps, props = zip(*fps_and_props)
        
        # Turn props into list of named tuples
        if props_named_tuple:
            props = [props_named_tuple(*p) for p in props]
        
        pickle.dump(props,propfile)
        pickle.dump(fps,fpfile)
        
    propfile.close()
    fpfile.close()
    
    print('Combining pickle files')
    sys.stdout.flush()
    props = []
    with open(os.path.join(outpath,"props.pickle"), 'rb') as pf:
        try:
            while True:
                props.extend(pickle.load(pf))
        except EOFError:
            pass
    with open(os.path.join(outpath,"props.pickle"), "wb+") as f:
        pickle.dump(props,f)
    
    fps = []
    with open(os.path.join(outpath,"fps.pickle"), 'rb') as pf:
        try:
            while True:
                fps.extend(pickle.load(pf))
        except EOFError:
            pass
    with open(os.path.join(outpath,"fps.pickle"), "wb+") as f:
        pickle.dump(fps,f)
    print('Batch processing complete!')
    sys.stdout.flush()
    return fps, props
    
def MolsToLSHForest(mol_list, save_path = "./", worker = os.cpu_count()-1, batch_size = None):
    
    print('Available CPU Cores =', os.cpu_count())
    print('Number of CPU Core used =', worker)
    print('\nTotal Number of Mols =', len(mol_list))
    if batch_size: print('Batch Size =', batch_size)
    if not os.path.exists(outpath): os.makedirs(outpath)
    print('Saving Files at', outpath)
    sys.stdout.flush()
    
    if batch_size:
        fps, props = batch_convert(mol_list, batch_size, outpath, props_named_tuple = Props)
    else:
        fps, props = single_convert(mol_list, outpath, props_named_tuple = Props)
            
    print("Loading data and converting to LSH Forest data")
    print('Converting MinHash Fingerprints to Vectors')
    sys.stdout.flush()
    fps = [tm.VectorUint(fp) for fp in fps]
    print(len(fps),'Fingerprints Converted')
    sys.stdout.flush()
    
    lf.batch_add(fps)
    lf.index()
    lf.store(os.path.join(outpath,"lf.dat"))
    
    return lf, props
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='LSH Forest Data Converter')
    parser.add_argument("filename", metavar = 'F', type=str, 
                        help="Path to the input file, can be either .smi (with no header) or .sdf")
    parser.add_argument("-b", "--batch", type=int, 
                        help="Size of the batch for processing, default to all (one batch)",
                       default = None)
    parser.add_argument("-o","--output", type=str, help="Output file path, default to ./LSHData", 
                        default = "./LSHData")
    parser.add_argument("-w","--worker", type=int, help="Number of workers (CPU cores) to use for multiprocessing,\
                        default to the number of available CPU cores minus one", 
                        default = os.cpu_count()-1)
    parser.add_argument("-d","--dim", type=int, help="Fingerprint dimension, default to 1024", 
                        default = 1024)
    
    a = parser.parse_args()
    outpath = os.path.abspath(a.output)
    mols = file_to_mols(a.filename)
    
    # Define a named properties tuple
    # To pickle a named tuple correctly:
    ## 1) The named tupple object has to be declared under __main__ 
    ## 2) The declared variable for the named tuple has to match 
    ##    the tuple name in the quotation mark!! 
    Props = namedtuple('Props',['SMILES','MolWt','LogP','QED','SAS'])
    
    # MinHash fingerprints (mhfp) encoder. This is a specialized molecular fingerprint scheme
    enc = MHFPEncoder(a.dim)
    # Locality Sensitive Hashing Forest
    lf = tm.LSHForest(a.dim, 64)
    
    MolsToLSHForest(mol_list = mols, 
                    save_path = outpath, 
                    worker = a.worker, 
                    batch_size = a.batch)
    