#!/usr/bin/env python
# -*- coding: utf-8 -*-

__doc__ = """Module for handling initial data formats and transforming them to requested by GraphMB tool:
 - contact map transformation to contig contact graph
 - abundance estimation from contig names (if available)
 - labels formatting for AMBER input convention
 - visuzualization of contig lengths distribution
  
 """

import os
import sys
import shutil

package_dir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, os.path.join(package_dir, "working_dir"))  # Add directory to path if User hasn`t done it yet

import logging
import argparse
import contact_map_processing
import io_prep_tools
from tqdm import tqdm
import numpy as np


def initialize_logger():
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s")

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)

    file_handler = logging.FileHandler(os.path.join(args.outdir, "preprocessing.log"), mode="w", encoding="utf-8")
    file_handler.setFormatter(formatter)
    file_handler.setLevel(logging.DEBUG)

    root.addHandler(console_handler)
    root.addHandler(file_handler)
    return root


def check_header(file):
    with open(file) as labels_file:
        line_1, line_2, line_3 = (labels_file.readline().strip() for _ in range(3))
        third_line = set(line_3.split("\t"))
        if any((line_1.split(":")[0] != "@Version",
                line_2.split(":")[0] != "@SampleID",
                not all(("@@SEQUENCEID" in third_line,
                         "BINID" in third_line,
                         "_LENGTH" in third_line))
                )):
            return False
        return True


if __name__ == '__main__':

    # Parse arguments:
    parser = argparse.ArgumentParser(description="Preprocessing Hi-C contact map for GraphMB or DMoN, \
    abundances for GraphMB, labels for AMBER")

    # Arguments for contact map
    parser.add_argument("-c", "--contact_map", type=str, default="contact_map.tsv", help="File with contact map")
    parser.add_argument("--score_column", type=str, default="SpadesScore",
                        help="Column with measure of Hi-C adjacency")
    # parser.add_argument("-d", "--")
    parser.add_argument("--scaling", type=str, choices=["sqrt", "log"], default="log",
                        help="Method of scaling Hi-C adjacency score")
    parser.add_argument("-e", "--edge_columns", type=str, nargs=2, default=["FirstName", "SecondName"],
                        help="Column names from contact map containing contig names with Hi-C link")
    parser.add_argument("--feature_columns", type=str, nargs="+", required=False,
                        help="""Features from contact map for DMoN (will be ignored if you use GraphMB)
                                Order must be: Feature_i_node_left, Feature_i_node_right, Feature_j_node_left, Feature_j_node_right etc.""")

    parser.add_argument("--graphmb", action="store_true", help="Process for GraphMB tools")
    parser.add_argument("--dmon", action="store_true", help="Process for DMoN tool")
    parser.add_argument("--vamb", action="store_true", help="Process for DMoN tool")

    # arguments for contigs
    parser.add_argument("-f", "--fasta_assembly", "--fasta", type=str, default="scaffolds.fasta",
                        help="File with assembled contigs/scaffolds")
    parser.add_argument("-a", "--abundances", type=str,
                        help="File with contig abundances -> transform for VAMB OR GraphMB")
    parser.add_argument("--mimic-jgi", "--mimic_jgi", action="store_true",
                        help="Estimate contigs coverage from contig names")

    parser.add_argument("-o", "--outdir", type=str,
                        default=os.path.join(os.getcwd(), "input"),
                        help="Directory for output to be written. MUST be empty. Or use --force to overwrite")
    parser.add_argument("--force", action="store_true",
                        help="Force override output dir if it isn`t empty. BE CAREFUL!")

    # arguments for labels
    parser.add_argument("--labels", type=str, help="File with known ground-truth")
    parser.add_argument("--header", type=str, help="Path for header file for AMBER tool")
    parser.add_argument("--draw", action="store_true", help="Option for drowing contigs lengths distribution")

    args = parser.parse_args()

    # Initializing logger

    root = initialize_logger()

    # Initiate output writing:

    if os.path.isdir(args.outdir) and not args.force:
        raise FileExistsError(
            f"Directory {args.outdir} already exists! Specify another one or use '--force' to overwrite")

    elif os.path.isdir(args.outdir) and args.force:
        shutil.rmtree(args.outdir, ignore_errors=True)
        # Inform user about directory overwriting:
        root.warning(f"Force overwriting directory {args.outdir}")

    os.mkdir(args.outdir)

    root.info(f"Reading contact map from {args.contact_map}")

    contact_map = contact_map_processing.ContactMap(
        path=args.contact_map,
        feature_columns=args.feature_columns,
        node_1_column=args.edge_columns[0],
        node_2_column=args.edge_columns[1],
        score_column=args.score_column,
        scaling_method=args.scaling,
    )

    if args.labels and not check_header(args.labels):
        root.info("Labels are incorrect for AMBER, loading from headers file...")
        if not check_header(args.header):
            msg = "Provided header file doesn't contain valid header for AMBER!"
            root.critical(msg)
            raise AssertionError(msg)

        labels_AMBER_dataset = io_prep_tools.AmberDataPreprocessor(
            data=args.labels,
            header_file=args.header,
        )
        output_path = os.path.join(args.outdir, "ground_truth.tsv")

        root.info(f"Saving labels with correct AMBER format in {output_path}")
        labels_AMBER_dataset.create_output_format(output_path)
        root.info("Done!")
    else:
        root.info("Provided labels are correctly organised for AMBER work, well done!")

    root.info("Observing contigs compound")
    scaffolds_full = []
    n_lines = sum(1 for i in open(args.fasta_assembly, 'rb'))
    with open(args.fasta_assembly) as scaffolds:
        for line in tqdm(scaffolds, total=n_lines):
            if line.startswith(">"):
                scaffolds_full.append(line[line.find(">") + 1:].strip())

    if args.graphmb:
        root.info("Started preparation for GraphMB requirements")
        root.info(f"Found {len(scaffolds_full)} contigs in {args.fasta_assembly}")

        contig2sequence = io_prep_tools.create_gfa_file(
            hic_data=contact_map,
            required_contigs=scaffolds_full,
            scaffolds_file=args.fasta_assembly,
            saving_dir=args.outdir,
            mimic=args.mimic_jgi and not args.abundances,
        )

        if args.draw and args.fasta_assembly:
            hic_scaffolds = contact_map.contigs_with_hic
            lengths_of_not_hic_scaffolds = np.array(
                [len(contig2sequence[scaffold]) for scaffold in scaffolds_full if scaffold not in hic_scaffolds])
            all_lengths = np.array(list(map(len, contig2sequence.values())))
            io_prep_tools.plot_contigs_lengths_distribution(all_lengths=all_lengths,
                                                            lengths_of_hic_scaffolds=lengths_of_not_hic_scaffolds,
                                                            outdir=args.outdir,
                                                            )

    if args.abundances:
        io_prep_tools.abundance2jgi(abundances_file=args.abundances,
                                    save_dir=args.outdir
                                    )

    if args.dmon:
        root.info("Creating adjacency & feature matrices for DMoN")
        sparse_adj, sparse_features = contact_map.get_sparse_adjacency_feature_matrices()

        root.info(f"Saving contact graph to {os.path.join(args.outdir, 'contact_map')}.npz")
        io_prep_tools.create_npz(out_file_name=os.path.join(args.outdir, "contact_map"),
                                 features_csr=sparse_features,
                                 adj_csr=sparse_adj)

    root.info(f"DONE! All results are saved in {args.outdir}")
