# Make TMAP with minhash fingerprint from the generated designs
# Plot with Faerun
# This is part of SampleDock package
#
# COMMERCIAL LICENSE: https://github.com/atfrank/SampleDock/blob/master/LICENSE_COMMERICAL
# NON-COMMERCIAL LICENSE: https://github.com/atfrank/SampleDock/blob/master/LICENSE_NON-COMMERICAL
# Frank Lab 2020-2021

import pickle
import tmap as tm
from mhfp.encoder import MHFPEncoder
from faerun import Faerun
from matplotlib import pyplot as plt
import pandas as pd

from tqdm.contrib.concurrent import process_map

def df_from_molProps(mols):
    # declare a named tuple
    Prop = namedtuple('Prop',['SMILES','Score','MolWeight','LogP','QED','SAS'])
    props = [(mol.GetProp('SMILES'),
              float(mol.GetProp('SCORE.INTER')),
              float(mol.GetProp('MolWeight')),
              float(mol.GetProp('LogP')),
              float(mol.GetProp('QED')),
              float(mol.GetProp('SAS'))) for mol in mols]
    # Make it a named tuple
    props = [Prop._make(p) for p in props]
    return pd.DataFrame(props)
        
def LSH_Convert(mols, num_workers):
    # MinHash fingerprints (mhfp) encoder for molecular fingerprinting
    enc = MHFPEncoder(1024)
    # Locality Sensitive Hashing Forest Instance
    lf = tm.LSHForest(1024, 64)
    
    print("Number of mols to be hashed:", len(mols))
    fps = process_map(enc.encode_mol, 
                      mols, 
                      chunksize = 100, 
                      max_workers=num_workers)
    
    fp_vecs = [tm.VectorUint(fp) for fp in fps]
    
    lf.batch_add(fp_vecs)
    lf.index()
    # save fp and lf   
    with open(os.path.join(outpath,"fps.pickle"), "wb") as fpfile:
        pickle.dump(fps,fpfile)
    lf.store(os.path.join(outpath,"lf.dat"))
    print('LSH data files saved!')
    return lf

def tree_coords(lf, node_size = 1/20, k = 20, mmm_rps = 2):
    print('Converting to tmap coordinates')
    # Create a LayoutConfiguration instance
    cfg = tm.LayoutConfiguration()
    cfg.node_size = node_size
    cfg.mmm_repeats = mmm_rps
    cfg.sl_extra_scaling_steps = 5
    cfg.k = k
    cfg.sl_scaling_type = tm.RelativeToAvgLength

    #Create minimum spanning tree from the LSHForest and LayoutConfiguration instance
    #The x and y coordinates of the vertices, the ids of the vertices spanning the edges
    #information on the graph is ignored

    x, y, s, t, _ = tm.layout_from_lsh_forest(lf, cfg)
    return list(x), list(y), list(s), list(t)

def df_to_faerun(df):
    print('Making Faerun plot')
    f = Faerun(view="front", coords=False)
    f.add_scatter(
        # No space in the string allowed for the name, use underscore!!
        # Cannot start with a number, it has to be a letter!! Weird Bug!!
        # My guess is that the string is to be converted to a variable name, 
        # therefore it has to be compatible with python variable naming scheme
        "SampleDock",
        {
            "x": x,
            "y": y,
            "c": [
                df['Score'],
                df['MolWeight'],
                df['LogP'],
                df['QED'],
                df['SAS']
            ],
            "labels": df['SMILES'],
        },
        shader="smoothCircle",
        point_scale=2.0,
        max_point_size=20,
        categorical=[False, False, False, False, False],
        colormap=["rainbow_r","rainbow", "rainbow", "rainbow", "Blues"],
        series_title=[
            "Docking Score",
            "Molecular Weight",
            "Lipophilicity",
            "Quantitative Estimate of Druglikeness",
            "Synthetic Accessibility Score",
        ],
        has_legend=True,
        )
    # The first character of the name has to be a letter!
    f.add_tree("SnD_Tree", {"from": s, "to": t}, point_helper="SampleDock")
    
    print('Plotting finished')
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='Plotter for tmap with Faerun frontend.\
                        Use LSH_batch.py to convert large dataset with limited memory,\
                        or convert .sdf files that need to be calculated for physicochemical properties.')
    parser.add_argument("mol_filename", metavar = 'F', type=str, 
                        help="Path to the the .sdf file of the generated designs")
    parser.add_argument("-o","--output", type=str, help="Output file path, default to ./LSHData", 
                        default = "./LSHData")
    parser.add_argument("-w","--worker", type=int, help="Number of workers (CPU cores) to use for multiprocessing,\
                        default to the number of available CPU threads minus one", 
                        default = os.cpu_count()-1)
    
    a = parser.parse_args()
    outpath = os.path.abspath(a.output)
    
    mols = [m for m in SDMolSupplier(mol_filepath) if m]
    df = df_from_molProps(mols)
    lf = LSH_Convert(mols, outpath)
    mols.clear()
    
    x, y, s, t = tree_coords(lf)
    df['x'] = x
    df['y'] = y
    df['s'] = s
    df['t'] = t
    df.to_csv(os.path.join(outpath,"props.csv"),index = False)
    df_to_faerun(df)