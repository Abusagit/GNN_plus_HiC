import argparse
from pathlib import Path
from tqdm import tqdm
import pandas as pd

from GNN_plus_HiC.graphmb_workflow.construct_graphmb_wf_config import create_json

def get_parser():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-o", "--outdir", help="Directory for contigs")
    parser.add_argument("-i")
    parser.add_argument("--delim")
    
    return parser



def generate(outdir, input_table, delimiter):
    
    config_table = pd.read_csv(input_table, index_col=0, delimiter=delimiter)
    
    ...
    
    
    
    
    
    