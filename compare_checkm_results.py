__doc__ = """"""

import plotly.graph_objects as go
import ast
import pandas as pd
import argparse
import os
from collections import defaultdict
from tqdm import tqdm


PLOT_NAME = "completeness_contamination_scatterplot.png"


class CheckMResult:
    def __init__(self, path, delim="\t"):
        self.path = path
        self.data = self.read_checkm(path, delim)

    @staticmethod
    def read_checkm(path, delim):
        bin_ids = []
        columns = defaultdict(list)

        numlines = sum(1 for _ in open(path, 'rb'))
        with open(path) as f_read:

            for line in tqdm(f_read, total=numlines):
                bin_id, values = line.split(delim)

                bin_ids.append(bin_id)

                values_dict = ast.literal_eval(values)

                for feature, value in values_dict.items():
                    columns[feature].append(value)

        dataframe = pd.DataFrame(index=bin_ids, data=columns)[
            ["marker lineage", "Completeness", "Contamination", "GC", "Genome size",
             "# scaffolds", "# contigs", "N50 (scaffolds)", "N50 (contigs)", "Mean scaffold length",
             "Mean contig length"]]

        dataframe.index.name = "binid"

        dataframe[["Completeness", "Contamination"]] /= 100

        dataframe["Purity"] = 1 / (1 + dataframe["Contamination"])

        dataframe["F1"] = 2 * dataframe["Completeness"] * dataframe["Purity"] / (
                dataframe["Completeness"] + dataframe["Purity"])

        return dataframe

    def get_HQ_genomes(self, completeness_lower_bound=0.95, purity_lower_bound=0.95):

        return self.data.query(
            f"Completeness > {completeness_lower_bound} & Purity > {purity_lower_bound}")


def plot_different_tools_results(tools_results_paths, tool_names, min_completeness, min_purity, outdir):
    fig = go.Figure()

    for result, tool_name in zip(tools_results_paths, tool_names):

        df = CheckMResult(path=result).get_HQ_genomes(completeness_lower_bound=min_completeness,
                                                      purity_lower_bound=min_purity)

        fig.add_trace(go.Scatter(y=df["Completeness"], x=df["Purity"],
                                 mode="markers",
                                 name=tool_name,
                                 marker_size=12,
                                 ))

    fig.update_traces(mode='markers', marker_line_width=1, marker_size=25
                      )

    fig.update_layout(yaxis_title="Completeness", xaxis_title="Purity",
                      font=dict(
                          family="Proxima Nova",
                          size=20,
                          color="Dark Blue"
                      ),
                      title="Completeness & Purity for different tools",
                      legend=dict(
                          yanchor="bottom",
                          y=0.01,
                          xanchor="left",
                          x=0.01, orientation="h"))

    fig.write_image(os.path.join(outdir, PLOT_NAME))


if __name__ == '__main__':
    parser = argparse.ArgumentParser("Plot tools Completeness/Purity scatterplot for chosen tools")

    parser.add_argument("-i", "--input_checkm", type=str, nargs="+",
                        help="Checkm results stored at <CheckM_out/storage/bin_stats_ext.tsv>")
    parser.add_argument("-l", "--labels", type=str, help="Corresponding labels for plot", nargs="+")

    parser.add_argument("-c", "--min-completeness", "--min_completeness", type=float, default=0.7,
                        help="Minimum Completeness of a bin for being taken into analysis")
    parser.add_argument("-p", "--min-purity", "--min_purity", type=float, default=0.7,
                        help="Minumal Purity of a bin for being taken into analysis")

    parser.add_argument("-o", "--outdir", type=str, help="Name of output file", default=os.getcwd())

    args = parser.parse_args()

    print(f"Arguments are: {args}")
    plot_different_tools_results(tools_results_paths=args.input_checkm,
                                 tool_names=args.labels,
                                 min_completeness=args.min_completeness,
                                 min_purity=args.min_purity,
                                 outdir=args.outdir)
    print(f"Result have been saved to {os.path.join(args.outdir, PLOT_NAME)}")
