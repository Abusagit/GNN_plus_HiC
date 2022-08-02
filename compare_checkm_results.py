#!/usr/bin/env python
# -*- coding: utf-8 -*-

__doc__ = """Informative comparative plots with binning metrics from several CheckM"""

import plotly.graph_objects as go
import ast
import pandas as pd
import argparse
import os
from collections import defaultdict
from tqdm import tqdm
from io_prep_tools import CheckMResult
from pathlib import Path

PLOT_NAME = "summary_{}_mindepth_{}.png"
HTML_NAME = "summary_{}_mindepth_{}.html"
TSV_NAME = "summary_mindepth_{}.tsv"
CSV_NAME = "summary_mindepth_{}.csv"


def plot_different_tools_results(tools_results_paths, tool_names, min_completeness, min_purity, outdir,
                                 depths,
                                 binnings,
                                 mindepth,
                                 ):
    def get_contigs_summary(data, binning_file):
        nonlocal contig_depths_lengths

        binning = pd.read_csv(binning_file, sep="\t", comment="@", index_col=0, header=None, dtype=str)

        HQ_bins = set(data.index.values)

        overall_len = overall_binned_contigs = hq_total_len = hq_total_num = hq_length_binned_in_hq_bin = hq_binned_in_hq_bin_num = 0

        requested_contigs = contig_depths_lengths.loc[binning.index, ["Length", "Depth"]].astype(int)
        binning["Length"], binning["Depth"] = requested_contigs["Length"], requested_contigs["Depth"]
        # for contig, (binid, *_) in tqdm(binning.iterrows(), total=binning.shape[0], desc="Checking binned contigs..."):
        #     # breakpoint()
        #     contig_length, contig_depth = contig_depths_lengths.loc[contig]
        #     overall_len += contig_length
            # contig_is_high_covered = contig_depth >= mindepth
            # contig_is_high_covered_and_in_hq_bin = contig_is_high_covered and binid in HQ_bins
            #
            # hq_total_num += contig_is_high_covered
            # hq_total_len += contig_is_high_covered * contig_length
            #
            # hq_binned_in_hq_bin_num += contig_is_high_covered_and_in_hq_bin
            # hq_length_binned_in_hq_bin += contig_is_high_covered_and_in_hq_bin * contig_length
            #
            # if contig_depth >= mindepth:
            #     hq_total_num += 1
            #     hq_total_len += contig_length
            #
            #     if binid in HQ_bins:
            #         hq_length_binned_in_hq_bin += contig_length
            #         hq_binned_in_hq_bin_num += 1

        overall_len = binning["Length"].sum()
        overall_binned_contigs = binning.shape[0]

        query = binning.query(f"Depth >= {mindepth}")
        hq_total_len = query["Length"].sum()
        hq_total_num = query.shape[0]

        filtered_index = list(filter(lambda c: binning.loc[c, 1] in HQ_bins, query.index.values))
        query = query.loc[filtered_index]
        hq_binned_in_hq_bin_num = query.shape[0]
        hq_length_binned_in_hq_bin = query["Length"].sum()

        # breakpoint()
        return (overall_len, hq_total_len, hq_length_binned_in_hq_bin,
                overall_binned_contigs, hq_total_num, hq_binned_in_hq_bin_num)

    contig_depths_lengths = pd.read_csv(depths, index_col=0,
                                        sep="\t").iloc[:, :2].rename(columns={"contigLen": "Length",
                                                                              "totalAvgDepth": "Depth"})

    fig = go.Figure()

    summary = []

    for i, (result, tool_name, binning_result) in enumerate(zip(tools_results_paths, tool_names, binnings)):
        df = CheckMResult(path=result).get_HQ_bins(completeness_lower_bound=min_completeness,
                                                   purity_lower_bound=min_purity
                                                   )

        summary.append([tool_name, df.shape[0],
                        *get_contigs_summary(data=df, binning_file=binning_result)])  # name, # of HQ genomes
        # breakpoint()

        fig.add_trace(go.Scatter(y=df["Completeness"], x=df["Purity"],
                                 mode="markers",
                                 name=f"{tool_name}\n({df.shape[0]})",
                                 marker_size=df["Genome size"] / 100_000,
                                 opacity=0.99,
                                 hovertext=df["Genome size"]

                                 ))

    summary = pd.DataFrame(data=summary,
                           columns=["Tool",
                                    "Total bins",
                                    "Overall binned length",
                                    "Length of binned to HQ bins",
                                    "Length of high-covered binned to HQ bins",
                                    "Overall # binned contigs",
                                    "Total # high-covered binned contigs",
                                    "# conigs binned to HQ bins",
                                    ])

    fig.update_traces(mode='markers', marker_line_width=1,
                      # marker_size=20,
                      # marker_opacity=0.9
                      )

    fig.update_layout(yaxis_title="Completeness", xaxis_title="Purity",
                      font=dict(
                          family="Proxima Nova",
                          size=20,
                          color="Dark Blue"
                      ),
                      title="Completeness & Purity of different tools",
                      legend=dict(
                          yanchor="top",
                          y=0.99,
                          xanchor="left",
                          x=0.01, orientation="h"))

    names = '_'.join(tool_names)

    fig.write_image(os.path.join(outdir, PLOT_NAME.format(names, mindepth)))
    fig.write_html(os.path.join(outdir, HTML_NAME.format(names, mindepth)))
    summary.to_csv(os.path.join(outdir, TSV_NAME.format(mindepth)), sep="\t", index=False)
    summary.to_csv(os.path.join(outdir, CSV_NAME.format(mindepth)), index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser("Plot tools Completeness/Purity scatterplot for chosen tools")

    parser.add_argument("-i", "--input_checkm", type=str, nargs="+",
                        help="Checkm results stored at <CheckM_out/storage/bin_stats_ext.tsv>")
    parser.add_argument("-l", "--labels", type=str, help="Corresponding labels for plot", nargs="+")

    parser.add_argument("-c", "--min-completeness", "--min_completeness", type=float, default=0.95,
                        help="Minimum Completeness of a bin for being taken into analysis")
    parser.add_argument("-p", "--min-purity", "--min_purity", type=float, default=0.95,
                        help="Minumal Purity of a bin for being taken into analysis")

    parser.add_argument("-b", "--binnings", type=str, nargs="+",
                        help="Files containing clusterings")

    parser.add_argument("-d", "--depths", type=str, help="File containing contig depths and lengths")
    parser.add_argument("--mindepth", type=int, help="Minimal depth for high-covered contig", default=100)

    parser.add_argument("-o", "--outdir", type=str, help="Name of output file", default=os.getcwd())

    args = parser.parse_args()

    print(f"Arguments are: {args}")

    Path(args.outdir).mkdir(exist_ok=True, parents=True)

    plot_different_tools_results(tools_results_paths=args.input_checkm,
                                 tool_names=args.labels,
                                 min_completeness=args.min_completeness,
                                 min_purity=args.min_purity,
                                 outdir=args.outdir,
                                 binnings=args.binnings,
                                 depths=args.depths,
                                 mindepth=args.mindepth)
    print(f"Result have been saved to {args.outdir}")
