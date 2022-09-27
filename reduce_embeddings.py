__doc__ = """Reduces dimensions of a given embeddings file in format [X, N] to format [X, New_dim]
X - # of elements
N - # of initial embeddings dimensions
New_dim - # of mapped dimensions

"""
import os

DEFAULT_PROCESSES = min(os.cpu_count(), 8)
os.environ['OPENBLAS_NUM_THREADS'] = str(DEFAULT_PROCESSES)
# These MUST be set before importing numpy
# I know this is a shitty hack, see https://github.com/numpy/numpy/issues/11826
os.environ["MKL_NUM_THREADS"] = str(DEFAULT_PROCESSES)
os.environ["NUMEXPR_NUM_THREADS"] = str(DEFAULT_PROCESSES)
os.environ["OMP_NUM_THREADS"] = str(DEFAULT_PROCESSES)


import argparse
from tqdm import tqdm
from pathlib import Path
import pickle
import numpy as np
from sklearn.manifold import TSNE
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
        tsne = TSNE(n_components=ndims, verbose=5)
        e_reduced = tsne.fit_transform(e.values)
        
        yield pd.DataFrame(e_reduced, columns=columns, index=e.index)
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser("Mapping embeddings to reduced vector space")
    parser.add_argument("-i", "--input", nargs="+", help="embeddings in .tsv format")
    parser.add_argument("-o", "--outdir", default=None)
    parser.add_argument("-d", "--ndimensions", default=3)
    
    args = parser.parse_args()
    
        
    if args.outdir:
        _dir = Path(args.outdir)
    else:
        _dir = Path("reduced_embeddings")
    
    outdir = _dir.parent.absolute() / _dir
    
    outdir.mkdir(exist_ok=True, parents=True)
    
    for e_new, old_name in zip(get_reduced_embeddings(*args.input, ndims=args.ndimensions), args.input):        
        
        old_name = '_'.join(os.path.splitext(old_name)[0].split('/'))
        old_path_absolute = Path(old_name).parent.absolute() / old_name
        
        new_name_path = DEFAULT_NAME.format(old_path_absolute)
        
        #save embedding to initial dir and create symlink located in outdir
        e_new.to_csv(new_name_path, sep='\t', index=True)
        os.symlink(src=new_name_path, dst=outdir / DEFAULT_NAME.format(old_name))
