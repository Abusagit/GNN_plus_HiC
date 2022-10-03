import argparse
from pathlib import Path
import json


def get_parser():
    
    parser = argparse.ArgumentParser("Automatic cretion of config with given parameters` values")
    
    parser.add_argument("-f", "--fasta", help="location of .fasta file (assembly)")
    parser.add_argument("-g", "--graph", help="location of assembly (contact) graph for GraphMB")
    parser.add_argument("-c", "--contact_map", help="location of contact map in .tsv format for embeddings aggregation")
    parser.add_argument("-p", "--processors", help="# of processors being used")
    parser.add_argument("--binners_root", "--binner_output_root", help="Root for Binners (VAMB. GraphMB) output dir" )
    parser.add_argument("--checkm_root", "--checkm_output_root", help="Root for CHECKM output dir")
    parser.add_argument("--aggregation_mode", help="Aggregation mode: full/during/after/none")
    parser.add_argument("-k", "--aggregation_neighbors", help="k nearest neighbors for aggregation")
    parser.add_argument("-i", "--aggregation_iterations", help="Iterations number for single iteration step")
    parser.add_argument("-m", "--contig_minlen", help="Threshold length for contig")
    parser.add_argument("-b", "--bin_minlen", help="Threshold length for bin in the output dir")
    parser.add_argument("--postfix_name")
    
    parser.add_argument("-o", "--output")
    
    
    return parser



def main():
    
    parser = get_parser()
    args = parser.parse_args()
    print(args)
    
    