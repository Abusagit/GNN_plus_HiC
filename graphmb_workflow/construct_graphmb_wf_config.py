import argparse
from pathlib import Path
import json



def get_parser():
    
    parser = argparse.ArgumentParser("Automatic cretion of config with given parameters` values")
    
    io_options = parser.add_argument_group(title="IO options", description="Technical options")
    io_options.add_argument("-o", "--output")
    
    config_options = parser.add_argument_group(title="Config options", description="Config options themself")
    
    config_options.add_argument("-f", "--fasta", help="location of .fasta file (assembly)")
    config_options.add_argument("-g", "--assembly_graph", help="location of assembly (contact) graph for GraphMB")
    config_options.add_argument("-d", "--depths_file", help="Contigs` depths file for VAMB")
    config_options.add_argument("-c", "--contact_map", help="location of contact map in .tsv format for embeddings aggregation")
    
    config_options.add_argument("-p", "--cores", help="# of processors being used", default="32")
    config_options.add_argument("--binner_output_root", help="Root for Binners (VAMB,  GraphMB) output dir", default="outputs/")
    config_options.add_argument("--checkm_output_root", help="Root for CHECKM output dir", default="checkm_reults/")
    config_options.add_argument("--aggregation_mode", help="Aggregation mode: full/during/after/none", default="full")
    config_options.add_argument("-k", "--aggregation_neighbors", help="k nearest neighbors for aggregation", default="30")
    config_options.add_argument("-i", "--aggregation_iterations", help="Iterations number for single iteration step", default="1")
    config_options.add_argument("-m", "--contig_minlen", help="Threshold length for contig", default="2000")
    config_options.add_argument("-b", "--bin_minlen", help="Threshold length for bin in the output dir", default="500000")
    config_options.add_argument("-s", "--short_threshold", help="Threshold for contig length to be considered in embeddings aggregation", default="100000")
    config_options.add_argument("--postfix_name", default='')
        
    technical_args = {"output", "o"}
    return parser, technical_args


def create_json(outfile, arguments, technical_args):
    config_arguments = dict((arg, arguments[arg]) for arg in arguments if arg not in technical_args)
    
    if '/' in outfile:
        directory = '/'.join(outfile.split('/')[:-1])
        Path(directory).mkdir(exist_ok=True, parents=True)
    with open(outfile, "w") as h:
        json.dump(config_arguments, h, indent=4, sort_keys=True)

def main():
    
    parser, technical_args = get_parser()
    args = vars(parser.parse_args())
    

        
    create_json(args["output"], args, technical_args)

        
    
if __name__ == "__main__":
    main()
