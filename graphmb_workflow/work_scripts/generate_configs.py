import argparse
from pathlib import Path
from tqdm import tqdm
import pandas as pd
import os

from GNN_plus_HiC.graphmb_workflow.construct_graphmb_wf_config import create_json, main

def get_parser():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-o", "--outdir", help="Directory for contigs")
    parser.add_argument("-i", "--input", help="Table with configs informations")
    parser.add_argument("--delimiter", help="Delimiter for a file", default=',')
    
    return parser



def generate(outdir, input_table, delimiter):
    
    config_table = pd.read_csv(input_table, index_col=0, delimiter=delimiter)
    
    keys = list(config_table.index.values)
    
    config_dicts = (dict(zip(keys, config_table.iloc[:, i])) for i in range(config_table.shape[1]))
    technical_args = {"output"}
    
    for config_d in tqdm(config_dicts, total=config_table.shape[1], desc="Creating configs from table"):
        outfile = os.path.join(outdir, config_d["output"])
        create_json(outfile=outfile, arguments=config_d, technical_args=technical_args)
        
        
def main():
    
    parser = get_parser()
    args = parser.parse_args()
    
    generate(outdir=args.outdir, input_table=args.input, delimiter=args.delimiter)
    

if __name__ == '__main__':
    main()