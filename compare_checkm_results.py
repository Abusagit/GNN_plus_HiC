#!/usr/bin/env python
# -*- coding: utf-8 -*-

__doc__ = """Informative comparative plots with binning metrics from several CheckM files"""

import plotly.graph_objects as go
from plotly.subplots import make_subplots
# import plotly.express as px

import pandas as pd
import argparse
import os
from collections import defaultdict
from tqdm import tqdm
from io_prep_tools import CheckMResult
from pathlib import Path

OVERALL_PLOT_NAME = "summary_{}_mindepth_{}.png"
OVERALL_HTML_NAME = "summary_{}_mindepth_{}.html"

TSV_NAME = "summary_{}_mindepth_{}.tsv"
CSV_NAME = "summary_{}_mindepth_{}.csv"

SUMMARY_PLOT_NAME = "columns_summary.{}"  # format [png, html]

FLOAT_FORMAT = "{:.2f} mbp".format


def plot_columns_from_summary(summary_df, tool_name_col, outdir):
    
    
    # summary_df.drop(tool_name_col, axis=1, inplace=True)
    summary_cols = summary_df.columns[1:] # tuple(filter(lambda x: x!= tool_name_col, summary_df.columns))

    traces = []
    titles = []
    # TODO make it prettier

    widths = [0.5 for _ in range(summary_df.shape[0])]
    for column in tqdm(summary_cols, total=summary_df.shape[1] - 1, desc="Drawing plots for summary"):
        title = column.replace('#', "Number")
        # breakpoint()
        trace = go.Bar(x=summary_df[tool_name_col], y=summary_df[column], width=widths, textposition='auto', 
        texttemplate="%{y:.2f}" if summary_df[column].dtype != "int64" else "%{y}",
        text=summary_df[column])
        traces.append(trace)
        titles.append(title)

    fig = make_subplots(rows=len(summary_cols), cols=1, shared_xaxes=False, vertical_spacing=0.08, 
    subplot_titles=titles)

    for subplot_index, trace in enumerate(traces, 1):
        fig.add_trace(trace, row=subplot_index, col=1)


    fig.update_layout(height=300 * len(summary_cols), width= 150 * (summary_df.shape[0]),
                  title_text="Tools different statistics",
                    # font=dict(
                    #       family="Proxima Nova",
                    #       size=20,
                    #       color="Dark Blue"
                    #   ),
                  showlegend=False)
    fig.write_html(os.path.join(outdir, SUMMARY_PLOT_NAME.format("html")))
    fig.write_image(os.path.join(outdir, SUMMARY_PLOT_NAME.format("png")))


def plot_different_tools_results(tools_results_paths, tool_names, min_completeness, min_purity, outdir,
                                 depths,
                                 binnings,
                                 mindepth,
                                 ):
    def get_contigs_summary(hq_data, binning_file):
        nonlocal contig_depths_lengths

        binning = pd.read_csv(binning_file, sep="\t", comment="@", index_col=0, header=None, dtype=str)
    

        HQ_bins = set(hq_data.index.values)

        # overall_len = overall_binned_contigs = hq_total_len = hq_total_num = hq_length_binned_in_hq_bin = hq_binned_in_hq_bin_num = 0
        try:
            requested_contigs = contig_depths_lengths.loc[binning.index, ["Length", "Depth"]].astype(float)
        except:
            breakpoint()
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
        hq_total_len = hq_data["Genome size"].sum()
        hq_total_num = query.shape[0]

        # breakpoint()
        filtered_index = list(filter(lambda c: binning.loc[c, 1] in HQ_bins, query.index.values))
        query = query.loc[filtered_index]
        hq_binned_in_hq_bin_num = query.shape[0]
        # breakpoint()
        hq_length_binned_in_hq_bin = query["Length"].sum()
        # breakpoint()

        # breakpoint()

        return (len(HQ_bins), overall_len, hq_total_len, hq_length_binned_in_hq_bin,
                overall_binned_contigs, hq_total_num, hq_binned_in_hq_bin_num)

    contig_depths_lengths = pd.read_csv(depths, index_col=0,
                                        sep="\t").iloc[:, :2].rename(columns={"contigLen": "Length",
                                                                              "totalAvgDepth": "Depth"})

    fig = go.Figure()

    summary = []

    for _, (result, tool_name, binning_result) in enumerate(zip(tools_results_paths, tool_names, binnings)):
        df = CheckMResult(path=result).get_HQ_bins(completeness_lower_bound=min_completeness,
                                                   purity_lower_bound=min_purity
                                                   )

        summary.append([tool_name, *get_contigs_summary(hq_data=df, binning_file=binning_result)])  # name, # of HQ genomes
        # breakpoint()

        df.reset_index(inplace=True)
        fig.add_trace(go.Scatter(y=df["Completeness"], x=df["Purity"],
                                 mode="markers",
                                 name=f"{tool_name}\n({df.shape[0]})",
                                 marker_size=df["Genome size"] / 100_000,
                                #  opacity=0.99,
                                 hovertext=df["binid"]

                                 ))

    summary = pd.DataFrame(data=summary,
                           columns=["Tool",
                                    "Total HQ bins",
                                    "Total length (mbp)",
                                    "HQ bins length (mbp)",
                                    f"Length of high-covered (> {mindepth}) binned to HQ bins (mbp)",
                                    "Sequences binned",
                                    "High-covered binned",
                                    "High-covered binned to HQ bins",
                                    ]).sort_values(by="Total HQ bins", ascending=False)


    summary.iloc[:, 2:5] /= 1_000_000
    
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

    

    plot_columns_from_summary(summary_df=summary, tool_name_col="Tool", outdir=outdir)
    fig.write_image(os.path.join(outdir, OVERALL_PLOT_NAME.format(names, mindepth)))
    fig.write_html(os.path.join(outdir, OVERALL_HTML_NAME.format(names, mindepth)))
    summary.to_csv(Path(outdir, TSV_NAME.format(names, mindepth)), sep="\t", index=False, float_format=FLOAT_FORMAT)
    summary.to_csv(Path(outdir, CSV_NAME.format(names, mindepth)), index=False, float_format=FLOAT_FORMAT)

def get_parser():
    parser = argparse.ArgumentParser("Plot tools Completeness/Purity scatterplot for chosen tools")
    
    parser.add_argument("-i", "--input_checkm", type=str, nargs="+",
                        help="Checkm results stored at <CheckM_out/.. (provide only main main checkm directories including storage/bin_stats_ext.tsv)")
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



    return parser

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    print(f"Arguments are: {args}")

    Path(args.outdir).mkdir(exist_ok=True, parents=True)

    
    checkm_inputs = [checkm_input if checkm_input.split('.')[-1] == "tsv" else f"{checkm_input}/storage/bin_stats_ext.tsv" for checkm_input in args.input_checkm]
    
    # tools_results = [Path(checkm_result, "storage/bin_stats_ext.tsv") for checkm_result in args.input_checkm]
    plot_different_tools_results(tools_results_paths=checkm_inputs,
                                 tool_names=args.labels,
                                 min_completeness=args.min_completeness,
                                 min_purity=args.min_purity,
                                 outdir=args.outdir,
                                 binnings=args.binnings,
                                 depths=args.depths,
                                 mindepth=args.mindepth)
    print(f"Result have been saved to {args.outdir}")
