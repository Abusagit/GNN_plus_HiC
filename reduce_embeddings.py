__doc__ = """Reduces dimensions of a given embeddings file in format [X, N] to format [X, New_dim]
X - # of elements
N - # of initial embeddings dimensions
New_dim - # of mapped dimensions

"""

import argparse
from tqdm import tqdm
from pathlib import Path
import pickle
import os
import pandas as pd

DEFAULT_NAME = "{}_reduced.tsv"


def handle_filename(filename):
    if filename.split('.')[-1] == "pkl" or filename.split('.')[-1] == "pickle": # pickle with contig names as key and values as embeddings
        with open(filename, "rb") as handler:
            _data = pickle.load(handler)
            data = pd.DataFrame.from_dict(_data, orient="index")
            
    elif filename.split('.')[-1] == "tsv":
        data = pd.read_csv(filename, header=None, index_col=0, sep='\t')
    
    elif filename.split('.')[-1] == "npz": #numpy compressed format with "arr_0" as key
        with np.load(filename) as handler:
            data = handler["arr_0"]
            
    return data


def get_reduced_embeddings(*embs_names, ndims=2) -> pd.DataFrame:
    
    if ndims == 2:
        columns = ["x", "y"]
    elif ndims == 3:
        columns = ["x", "y", "z"]
    else:
        columns = list(range(ndims))
    
    
    for e_name in tqdm(embs_names, desc=f"Reducing dimensions to {ndims} for every given embedding file..."):
        e = handle_filename(e_name)
        tsne = TSNE(n_components=ndims, verbose=5, n_jobs=NJOBS)
        e_reduced = tsne.fit_transform(e.values)
        
        yield pd.DataFrame(e_reduced, columns=columns, index=e.index)
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser("Mapping embeddings to reduced vector space")
    parser.add_argument("-i", "--input", nargs="+", help="embeddings in .tsv format")
    parser.add_argument("-o", "--output", default=None)
    parser.add_argument("-d", "--ndimensions", default=3)
    parser.add_argument("--njobs", default=8, type=int)
    
    args = parser.parse_args()
    
    NJOBS = args.njobs
    
    names = [name.split('/')[-1].split('.')[0] for name in args.input]
    
    outdir = Path(*args.input[0].split("/")[:-1]) if '/' in args.input[0] else Path().cwd()
    
    outdir.mkdir(exist_ok=True, parents=True)
    
    os.environ['OPENBLAS_NUM_THREADS'] = str(args.njobs)
    
    import numpy as np
    from sklearn.manifold import TSNE


    if not args.output:
        
        tmp = Path("reduced_embeddings")
        tmp.mkdir(exist_ok=True, parents=True)
        
        for e_new, e_old in zip(get_reduced_embeddings(*args.input, ndims=args.ndimensions), names):
            outname = str(outdir / DEFAULT_NAME.format(e_old))
            e_new.to_csv(outname, sep='\t', index=True)
            
            os.symlink(tmp / outname, tmp / DEFAULT_NAME.format(e_old))
            
    else:
        for e_new, e_old in zip(get_reduced_embeddings(*args.input, ndims=args.ndimensions), names):
            outname = str(outdir / DEFAULT_NAME.format(e_old))
            e_new.to_csv(outname, sep='\t', index=True)