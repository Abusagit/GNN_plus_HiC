#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import plotly.express as px
from tqdm import tqdm
import logging
from collections import defaultdict
import os
import sys
import pandas as pd
import numpy as np
import re

package_dir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, os.path.join(package_dir, "working_dir"))

import contact_map_processing

if __name__ == '__main__':
    DEFAULT_OUTDIR = "./"

    parser = argparse.ArgumentParser(description="Draw heatmap of Hi-C connections between and within clusters")

    parser.add_argument("-b", "--binning_results", help="Path to binning results",
                        type=str)
    parser.add_argument("-g", "--ground_truth", type=str,
                        help="Path to golden standard")
    parser.add_argument("-c", "--contact_map", type=str,
                        help="Path to contact map")
    parser.add_argument("-a", "--amber_summary", "--amber-summary", type=str, help="AMBER summary file")

    # parser.add_argument("-d", "--depths", type=str, help="Depths of contigs")
    parser.add_argument("-o", "--outdir", default=DEFAULT_OUTDIR)

    args = parser.parse_args()

    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s")

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.DEBUG)
    console_handler.setFormatter(formatter)
    root.addHandler(console_handler)
    #
    # root.info(f"Reading depths file {args.depths}")
    # depths_ = pd.read_csv(args.depths, sep="\t").query("contigLen >= 2000")
    # contig_depths = {contig: depth for contig, depth in
    #                  zip(depths_["contigName"].values, depths_["totalAvgDepth"].values)}
    # del depths_

    root.info(f"Reading contact map from {args.contact_map}")
    contact_map = contact_map_processing.ContactMap(path=args.contact_map)

    root.info("Obtaining Hi-C connections")
    hic_connections = set()
    for c_1, c_2 in tqdm(contact_map.data[["FirstName", "SecondName"]].values,
                         total=contact_map.data.shape[0]):
        hic_connections.add((c_1, c_2))
        hic_connections.add((c_2, c_1))

    root.info(f"Loading binning results from {args.binning_results}")
    labels = pd.read_csv(args.binning_results, skiprows=2, sep="\t")

    root.info("Assigning each contig to cluster")
    contig2bin = {}
    for _, row in tqdm(labels.iterrows(), total=labels.shape[0]):
        contig2bin[row["@@SEQUENCEID"]] = row["BINID"]

    root.info(f"Assigning each contig to corresponding genome according to ground truth {args.ground_truth}")

    contig2genome = {}
    with open(args.ground_truth) as f_read:
        for line in f_read:
            if line.startswith("@"):
                continue
            contig, genome, *_ = line.strip().split("\t")

            contig2genome[contig] = genome

    if args.amber_summary:
        genome_lengths = {}
        bin2genome = {}
        root.info(f"Loading AMBER summary from {args.amber_summary}")
        amber_summary = pd.read_csv(args.amber_summary,
                                    sep="\t")

        gold_standard = amber_summary.query("Tool == 'Gold standard'")
        for _, row in gold_standard.iterrows():
            genome_lengths[row["genome_id"]] = row["tp_length"]

        amber_summary = amber_summary.query("Tool == 'GraphMB'").set_index("BINID")
        amber_summary.index = amber_summary.index.astype(int)
        for binid, row in amber_summary.iterrows():
            bin2genome[binid] = row["genome_id"]

        root.debug(amber_summary)

    root.info("Building data for heatmap")

    unvisited_connections = hic_connections.copy()
    clustered_contigs = set(contig2bin.keys())

    HiC_LINKS_BETWEEN_BINS = defaultdict(lambda: defaultdict(int))
    HIC_links_between_contigs = defaultdict(list)
    i = 0
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

        i += 1 + cluster_1 != cluster_2
        if i % 1000 == 0:
            root.info(f"Processed {i} pairs")
    #
    # root.info("Estimating high-covered contigs misbinnings...")
    #
    # percentage_of_genome_high_covered_misclustered_HIC_LINKS = defaultdict(float)
    # percentage_of_genome_high_covered_misclustered_NO_HIC_LINKS = defaultdict(float)
    #
    # for _, (contig, cluster, contig_length) in tqdm(labels.iterrows(), total=labels.shape[0]):
    #     contig_genome = contig2genome.get(contig)
    #     if not contig_genome:
    #         continue
    #     ratio_of_genome = contig_length / genome_lengths[contig_genome]
    #     contig_coverage = contig_depths.get(contig) if contig_depths.get(contig) else float(re.findall("cov_(\d+\.\d*)", contig)[0])
    #
    #     if contig_coverage < 1000 and ratio_of_genome < 0.01:
    #         continue
    #
    #     contig_binid = contig2bin[contig]
    #     try:
    #         genome_in_contig_bin = amber_summary.loc[contig_binid, "genome_id"]
    #         genome_length = genome_lengths[genome_in_contig_bin]
    #         ratio_of_genome_in_bin = amber_summary.loc[contig_binid, "tp_length"] / genome_length
    #
    #         bin_size_sequences = amber_summary.loc[contig_binid, "total_seq_counts"]
    #     except KeyError:
    #         continue
    #
    #     if bin_size_sequences > 1:
    #         continue
    #
    #     for contig_2 in HIC_links_between_contigs[contig]:
    #         if contig_2 not in contig2genome:
    #             continue
    #
    #         try:
    #             contig2_binid = contig2bin[contig_2]
    #             contig2_genome = contig2genome[contig_2]
    #             genome_in_contig2_bin = amber_summary.loc[contig2_binid, "genome_id"]
    #             contig2_bin_TP_length = amber_summary.loc[contig2_binid, "tp_length"] / genome_length
    #             contig2_bin_TP_size_sequences = amber_summary.loc[contig2_binid, "tp_seq_counts"]
    #
    #             if all((
    #                     contig_genome == contig2_genome, contig_binid != contig2_binid, ratio_of_genome_in_bin >= 0.5
    #             )):
    #                 percentage_of_genome_high_covered_misclustered_HIC_LINKS[contig_genome] += ratio_of_genome
    #                 break
    #         except KeyError:
    #             continue
    #     else:
    #         if HIC_links_between_contigs[contig] and not genome_length == contig_length:
    #             percentage_of_genome_high_covered_misclustered_NO_HIC_LINKS[contig_genome] += ratio_of_genome
    #
    # breakpoint()

    all_binids = sorted(filter(lambda x: all((sum(HiC_LINKS_BETWEEN_BINS[x].values()) - HiC_LINKS_BETWEEN_BINS[x][x] > 10,
                                              x in bin2genome)),
                               HiC_LINKS_BETWEEN_BINS.keys()), key=lambda bin_: bin2genome[bin_])

    root.info("Creating dataframe...")
    data_for_heatmap = np.array([
        [HiC_LINKS_BETWEEN_BINS[binid].get(bin_id_internal, 0) for bin_id_internal in all_binids] for binid in reversed(all_binids)])

    root.debug(f"{data_for_heatmap.shape=}")

    root.info("Drawing heatmap")

    index = set(amber_summary.index)
    all_binids = tuple(
        map(lambda x: f"{x} ({amber_summary.loc[x, 'genome_id'] if x in index else 'None'})", all_binids))

    fig = px.imshow(np.log10(data_for_heatmap + 1),
                    labels=dict(x="BINID", y="BINID", color="Log10(Hi-C links)"),
                    x=all_binids,
                    y=all_binids[::-1],
                    aspect="auto",
                    color_continuous_scale="BuPu",
                    height=1500,
                    width=1500,
                    )

    z_text = [
        [f"{data_for_heatmap[i, j]} Hi-C links" for j in range(data_for_heatmap.shape[1])] for i in range(
            data_for_heatmap.shape[0])
    ]

    # fig.update_traces(hovertemplate="Bins %{x} X %{y}<extra></extra>")
    # fig.update_xaxes(visible=False)
    # fig.update_yaxes(visible=False)
    fig.write_image(os.path.join(args.outdir, "heatmap.png"))
    fig.write_html(os.path.join(args.outdir, "heatmap.html"))

    root.info(f"Heatmap saved at {args.outdir}")
