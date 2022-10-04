import argparse
from pathlib import Path
import json


def get_parser():
    
    parser = argparse.ArgumentParser("Automatic cretion of config with given parameters` values")
    
    io_options = parser.add_argument_group(title="IO options", description="Technical options")
    io_options.add_argument("-o", "--output")
    
    config_options = parser.add_argument_group(title="Config options", description="Config options themself")
    
    config_options.add_argument("-f", "--fasta", help="location of .fasta file (assembly)")
    config_options.add_argument("-g", "--graph", help="location of assembly (contact) graph for GraphMB")
    config_options.add_argument("-c", "--contact_map", help="location of contact map in .tsv format for embeddings aggregation")
    config_options.add_argument("-p", "--processors", help="# of processors being used")
    config_options.add_argument("--binners_root", "--binner_output_root", help="Root for Binners (VAMB,  GraphMB) output dir" )
    config_options.add_argument("--checkm_root", "--checkm_output_root", help="Root for CHECKM output dir")
    config_options.add_argument("--aggregation_mode", help="Aggregation mode: full/during/after/none")
    config_options.add_argument("-k", "--aggregation_neighbors", help="k nearest neighbors for aggregation")
    config_options.add_argument("-i", "--aggregation_iterations", help="Iterations number for single iteration step")
    config_options.add_argument("-m", "--contig_minlen", help="Threshold length for contig")
    config_options.add_argument("-b", "--bin_minlen", help="Threshold length for bin in the output dir")
    config_options.add_argument("--postfix_name")
        
    
    return parser



def main():
    
    parser = get_parser()
    args = parser.parse_args()
    print(args)
    
    
if __name__ == "__main__":
    main()
