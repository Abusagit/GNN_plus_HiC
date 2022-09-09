#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from ast import parse
import plotly.express as px
from tqdm import tqdm
import logging
from collections import defaultdict
import os
import sys
import pandas as pd
import numpy as np
import re
from pathlib import Path

package_dir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, os.path.join(package_dir, "working_dir"))

import contact_map_processing
from io_prep_tools import CheckMResult


def get_parser():
    parser = argparse.ArgumentParser(description="Draw heatmap of Hi-C connections between and within clusters")

    parser.add_argument("-b", "--binning_results", help="Path to binning results",
                        type=str)
    parser.add_argument("-g", "--ground_truth", type=str,
                        help="Path to golden standard FOR AMBER ONLY")

    parser.add_argument("-m", "--contact_map", type=str,
                        help="Path to contact map")
    parser.add_argument("-s", "--scale_score",  type=str, help="Scaling function for Hi-C score",
                        choices=["sqrt", "log"], default=None)

    parser.add_argument("--tool", type=str, choices=["checkm", "amber"])
    parser.add_argument("-r", "--report", type=str, help="AMBER or CHECKM summary file")

    parser.add_argument("-d", "--depths", type=str, help="Depths of contigs")
    parser.add_argument("--score_norm", choices=["length", "num"], default=None)
    parser.add_argument("--lognorm_score", action="store_true")

    parser.add_argument("-o", "--outdir", default=os.getcwd())
    parser.add_argument("-l", "--label", help="Label of a plot, e.g.: <label>_heatmap.png")

    parser.add_argument("--only_hq", "--only-hq", action="store_true")
    parser.add_argument("--report_only_counts", action="store_true")

    return parser


def get_logger():
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s")

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.DEBUG)
    console_handler.setFormatter(formatter)
    root.addHandler(console_handler)

    return root


if __name__ == '__main__':

    parser = get_parser()
    args = parser.parse_args()
    root = get_logger()
    
    root.info(args)
    root.info(f"Reading contact map from {args.contact_map}")
    contact_map = contact_map_processing.ContactMap(path=args.contact_map, scaling_method=args.scale_score)

    root.info("Obtaining Hi-C connections")
    
    hic_pairs = set()
    hic_connections = defaultdict(lambda: defaultdict(float))
    
    for c_1, c_2, score in tqdm(contact_map.data[["FirstName", "SecondName", "SpadesScore"]].values,
                        total=contact_map.data.shape[0]):
        hic_connections[c_1][c_2] += score
        hic_connections[c_2][c_1] += score
        
        hic_pairs.add((c_1, c_2))
        hic_pairs.add((c_2, c_1))
        
    root.info(f"Loading binning results from {args.binning_results}")
    labels = pd.read_csv(args.binning_results, comment="@", sep="\t",
                        index_col=0, header=None, dtype=str)  # contigname | clusterr | [etc.]

    root.info("Assigning each contig to cluster")

    contig2bin = dict(zip(labels.index, labels.iloc[:, 0]))
    bin2seq_size = labels.groupby(1).size()
    bin2seq_size = dict(zip(bin2seq_size.index, bin2seq_size.values))

    if args.tool == "amber":
        bin2genome = {}
        root.info(f"Loading AMBER summary from {args.amber_summary}")
        amber_summary = pd.read_csv(args.amber_summary,
                                    sep="\t")

        gold_standard = amber_summary.query("Tool == 'Gold standard'")
        for _, row in gold_standard.iterrows():
            genome_lengths[row["genome_id"]] = row["tp_length"]

        amber_summary = amber_summary.query("Tool == 'GraphMB'").set_index("BINID")
        amber_summary.index = amber_summary.index  # .astype(int)
        for binid, row in amber_summary.iterrows():
            bin2genome[binid] = row["genome_id"]

        root.debug(amber_summary)
    else:  # CHECKM
        checkm_result = CheckMResult(args.report)
        if args.only_hq:
            hq_bins = checkm_result.get_HQ_bins()
            bin2genome = dict(zip(hq_bins.index, hq_bins["marker lineage"]))
            bin2length = dict(zip(hq_bins.index, hq_bins["Genome size"]))
        else:
            bin2genome = dict(zip(checkm_result.data.index, checkm_result.data["marker lineage"]))
            bin2length = dict(zip(checkm_result.data.index, checkm_result.data["Genome size"]))

    root.info(f"Assigning each contig to corresponding genome according to ground truth {args.ground_truth}")

    contig2genome = {}

    if args.tool == "amber":
        with open(args.ground_truth) as f_read:
            for line in f_read:
                if line.startswith("@"):
                    continue
                contig, genome, *_ = line.strip().split("\t")

                contig2genome[contig] = genome
    else:  # CHECKM
        for contig, bin_ in contig2bin.items():
            if bin_ in bin2genome:
                contig2genome[contig] = bin2genome[bin_]

    root.info("Building data for heatmap")

    unvisited_connections = hic_pairs.copy()
    clustered_contigs = set(contig2bin.keys())

    HiC_LINKS_BETWEEN_BINS = defaultdict(lambda: defaultdict(float))
    HIC_links_between_contigs = defaultdict(list)
    i = 0

    if args.report_only_counts:
        while unvisited_connections:
            c_1, c_2 = unvisited_connections.pop()
            unvisited_connections.discard((c_2, c_1))

            cluster_1 = contig2bin.get(c_1)
            cluster_2 = contig2bin.get(c_2)

            if not all((cluster_1, cluster_2)):
                continue
            HiC_LINKS_BETWEEN_BINS[cluster_1][cluster_2] += 1
            HiC_LINKS_BETWEEN_BINS[cluster_2][cluster_1] += cluster_1 != cluster_2

            HIC_links_between_contigs[c_1].append(c_2)
            HIC_links_between_contigs[c_2].append(c_1)

            i += 1 + (cluster_1 != cluster_2)
            if i % 1000 == 0:
                root.info(f"Processed {i} pairs")
    else:
        while unvisited_connections:
            c_1, c_2 = unvisited_connections.pop()
            unvisited_connections.discard((c_2, c_1))

            cluster_1 = contig2bin.get(c_1)
            cluster_2 = contig2bin.get(c_2)
            
            score = hic_connections[c_1][c_2]

            if not all((cluster_1, cluster_2)):
                continue
            HiC_LINKS_BETWEEN_BINS[cluster_1][cluster_2] += score
            HiC_LINKS_BETWEEN_BINS[cluster_2][cluster_1] += (cluster_1 != cluster_2) * score

            HIC_links_between_contigs[c_1].append(c_2)
            HIC_links_between_contigs[c_2].append(c_1)

            i += 1 + (cluster_1 != cluster_2)
            if i % 1000 == 0:
                root.info(f"Processed {i} pairs")



        all_binids = tuple(sorted(
            filter(lambda x: all((sum(HiC_LINKS_BETWEEN_BINS[x].values()) - HiC_LINKS_BETWEEN_BINS[x][x] > 10 * bool(args.report_only_counts),
                                x in bin2genome and bin2genome[x] != "root")),
        HiC_LINKS_BETWEEN_BINS.keys()), key=lambda bin_: bin2genome[bin_]))

    root.info("Creating dataframe...")
    
    if args.lognorm_score:
        raw_data = np.array([
            [np.log10(HiC_LINKS_BETWEEN_BINS[binid].get(bin_id_internal, 0) + 1) for bin_id_internal in all_binids]
            for binid in reversed(all_binids)]).astype(np.float64)
        
    else:
        raw_data = np.array([
            [HiC_LINKS_BETWEEN_BINS[binid].get(bin_id_internal,0) for bin_id_internal in all_binids]
            for binid in reversed(all_binids)]).astype(np.float64)
        
    data_for_heatmap = np.zeros_like(raw_data)
    
    if args.score_norm == 'length':
        for i in range(data_for_heatmap.shape[0]):
            binid_i = all_binids[len(all_binids) - i - 1]
            len_i = bin2length[binid_i]
            for j in range(i + 1, data_for_heatmap.shape[0]):
                binid_j = all_binids[j]
                len_j = bin2length[binid_j]
                
                data_for_heatmap[i, j] = raw_data[i, j] / ((binid_i + binid_j) // 100_000)
                data_for_heatmap[j, i] = raw_data[j, i] / ((binid_i + binid_j) // 100_000)
                
        for i in range(data_for_heatmap.shape[0]):
            len_i = bin2length[all_binids[len(all_binids) - i - 1]]
            data_for_heatmap[i, i] = raw_data[i,i] / (len_i // 100_000)
                
    
    elif args.score_norm == 'number':
        for i in range(data_for_heatmap.shape[0]):
            binid_i = all_binids[len(all_binids) - i - 1]
            seqs_num_i = bin2seq_size[binid_i]
            for j in range(i + 1, data_for_heatmap.shape[0]):
                binid_j = all_binids[j]
                seqs_num_j = bin2seq_size[binid_j]
                
                data_for_heatmap[i, j] = raw_data[i, j] / (seqs_num_i + seqs_num_j)
                data_for_heatmap[j, i] = raw_data[j, i] / (seqs_num_i + seqs_num_j)
                
        for i in range(data_for_heatmap.shape[0]):
            seqs_num_i = bin2seq_size[all_binids[len(all_binids) - i - 1]]
            data_for_heatmap[i, i] = raw_data[i,i] / seqs_num_i
    else:
        data_for_heatmap = raw_data
            
    root.debug(f"{data_for_heatmap.shape=}")

    root.info("Drawing heatmap")

    # if args.tool == "amber":
    #     index = set(amber_summary.index)
    # else:
    #     index = set(checkm_result.index)
    all_binids = tuple(
        map(lambda x: f"{x} ({bin2genome.get(x)}) LEN {bin2length[x]}", all_binids))
    
    if args.lognorm_score:
        COLOR = "Normalised Log10(Hi-C score)" if not args.report_only_counts else "Normalized Log10(Hi-C links)"
    else:
        COLOR = "HiC-score"
    fig = px.imshow(
        # np.log10(data_for_heatmap + 1),
        data_for_heatmap,
        labels=dict(x="BINID", y="BINID", color=COLOR),
        x=all_binids,
        y=all_binids[::-1],
        aspect="auto",
        color_continuous_scale="BuPu",
        height=1500,
        width=1500,
    )

    if args.report_only_counts:
        
        z_text = [
            [f"{data_for_heatmap[i, j]} normalised Hi-C links" for j in range(data_for_heatmap.shape[1])] for i in range(
                data_for_heatmap.shape[0])
        ]
    else:
        z_text = [
            [f"{data_for_heatmap[i, j]} normalised Hi-C score" for j in range(data_for_heatmap.shape[1])] for i in range(
                data_for_heatmap.shape[0])
        ]

    # fig.update_traces(hovertemplate="Bins %{x} X %{y}<extra></extra>")
    # fig.update_xaxes(visible=False)
    # fig.update_yaxes(visible=False)

    Path(args.outdir).mkdir(exist_ok=True, parents=True)
    fig.write_image(os.path.join(args.outdir, f"{args.label + '_' if args.label else ''}{'links_' if args.report_only_counts else 'scores_'}{'hq_' if args.only_hq else ''}heatmap.png"))
    fig.write_html(os.path.join(args.outdir, f"{args.label + '_' if args.label else ''}{'links_' if args.report_only_counts else 'scores_'}{'hq_' if args.only_hq else ''}heatmap.html"))

    root.info(f"Heatmap saved at {args.outdir}")
