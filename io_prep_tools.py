#!/usr/bin/env python
# -*- coding: utf-8 -*-

__doc__ = """Various functions for transforming input and output of programs"""

import logging
from tqdm import tqdm
import pandas as pd
import numpy as np
import os
import re
from functools import partial
from collections import defaultdict
import plotly.express as px
import bisect
import ast


logger = logging.getLogger(__name__)


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

    def get_HQ_bins(self, completeness_lower_bound=0.90, purity_lower_bound=0.90):

        return self.data.query(
            f"Completeness > {completeness_lower_bound} & Purity > {purity_lower_bound}")


class AmberDataPreprocessor:
    def __init__(self, data, header_file):
        self.header = self.store_header(header_file)
        self.data = pd.read_csv(data, sep="\t", index_col=0, header=0)

    @staticmethod
    def store_header(file):
        with open(file) as f_read:
            return ''.join(f_read.readline() for _ in range(3))

    def create_output_format(self, output_name):
        self.data.to_csv(f"{output_name}", sep="\t", index=True)

        with open(f"{output_name}", "r+") as f:
            content = f.read()
            f.seek(0, 0)
            f.write(self.header.rstrip("\r\n") + "\n" + content)


def create_npz(out_file_name, features_csr, adj_csr, encoder=None, labels=None, label_indices=None):
    npz = {"feature_data": features_csr.data,
           "feature_indices": features_csr.indices,
           "feature_indptr": features_csr.indptr,
           "feature_shape": np.array(features_csr.shape),
           "adj_data": adj_csr.data,
           "adj_indices": adj_csr.indices,
           "adj_indptr": adj_csr.indptr,
           "adj_shape": adj_csr.shape,
           "encoder_classes_": encoder.classes_ if encoder else np.arange(features_csr.shape[0]),
           "labels": labels,
           "label_indices": label_indices}

    np.savez(out_file_name, **npz)
    logger.info(f"Features file saved with name '{out_file_name}.npz'")


def transform_tnf_to_features(tnf_data: pd.DataFrame, encoder):
    nodes_from_contact_map = list(set(tnf_data.index) & set(encoder.classes_))
    nodes_encoded = encoder.transform(nodes_from_contact_map)

    shape = (len(encoder.classes_), tnf_data.shape[1])

    tnf_features = scipy.sparse.lil_matrix(np.zeros(shape))

    for node, node_index in zip(nodes_from_contact_map, nodes_encoded):
        node_tnf_data = tnf_data.loc[node]

        tnf_features[node_index] = node_tnf_data

    return scipy.sparse.csr_matrix(tnf_features)


def get_concatenated_features(old_features, additional_features):
    return scipy.sparse.hstack([old_features, additional_features])


def get_dict_from_npz(filepath):
    with np.load(filepath, "rb") as f_read:
        loader = dict(f_read)

    return loader


def mimic_jgi(contigs_with_info_in_names, initial_assembly_graph,
              # contact_map,
              lengths=None, save_name="features.txt",
              write_option=True):
    """
    Only suitable for data with the specidic pattern
    """

    len_cov_pattern = re.compile(r"length_(\d+)_cov_(\d+\.\d*)")

    regexp_vectorized = np.vectorize(lambda x: len_cov_pattern.findall(x)[0])

    class RegexpError(Exception):
        def __init__(self):
            self.message = "No information was found either in contig names or in assembly_graph. " \
                           f"Perhaps you need to add 'dp:i:' tag to {initial_assembly_graph}"
            super().__init__(self.message)

    try:
        if not lengths:
            lengths, coverages = regexp_vectorized(contigs_with_info_in_names)
        else:
            _, coverages = regexp_vectorized(contigs_with_info_in_names)
    except IndexError:  # re.findall has empty list <==> incorrect contig names
        logger.warning("Names of contigs don't provide lengths and depths, trying to fetch from contact map and gfa")
        try:
            logger.warning("Trying to compute coverage from edges")
            edge_coverage_pattern = re.compile("^S\t(\S+).*dp:i:(\d+)")  # (contig_name, depth)
            path_pattern = re.compile("^P\t(\S*)\t(\S+)")  # (contig, edges from path

            edge2coverage = {}
            contig2edges = defaultdict(list)

            with open(initial_assembly_graph) as gfa:
                for line in tqdm(gfa, desc="Reading GFA"):
                    edge_cov_p = edge_coverage_pattern.findall(line)
                    path_p = path_pattern.findall(line)
                    if edge_cov_p:
                        edge = edge_cov_p[0][0]
                        edge_coverage = int(edge_cov_p[0][1])
                        edge2coverage[edge] = edge_coverage
                        # breakpoint()
                    elif path_p:
                        contig = path_p[0][0]
                        path = list(map(lambda x: x[:-1], path_p[0][1].split(",")))
                        contig2edges[contig].extend(path)
                        # breakpoint()

            # if not lengths:
            #     lengths = [contact_map.contig2length[contig] for contig in contigs_with_info_in_names]

            coverages = []
            for index, contig in tqdm(enumerate(contigs_with_info_in_names), desc="Assgning depths to contigs"):
                contig_summarized_depth = sum(map(lambda e: edge2coverage[e], contig2edges[contig]))
                contig_length = lengths[index]

                coverages.append(contig_summarized_depth)
            breakpoint()
            assert len(coverages) == len(lengths) == len(contigs_with_info_in_names)

        except IndexError:
            breakpoint()
            raise RegexpError

    else:
        columns = ["contigName", "contigLen", "totalAvgDepth", "nodes.bam"]
        data = dict(zip(columns,
                        (contigs_with_info_in_names,
                         lengths,
                         coverages,
                         coverages,
                         )))

        dataframe = pd.DataFrame(data)
        if write_option:
            dataframe.to_csv(save_name, sep="\t", index=False)

            logger.info(f"Abundances were saved to {save_name}")

        return dataframe


def create_gfa_file(hic_data, required_contigs, scaffolds_file,
                    # contact_map,
                    initial_assembly_graph,
                    saving_dir="./", mimic=True,
                    save_gfa_file="assembly_graph.gfa",
                    save_fasta_file="assembly.fasta"):
    """    """

    # Processing scaffolds first:

    required_contigs_set = set(required_contigs)
    save_gfa_file = os.path.join(saving_dir, save_gfa_file)
    save_fasta_file = os.path.join(saving_dir, save_fasta_file)

    with open(save_gfa_file, "w") as save_gfa, open(save_fasta_file, "w") as save_fasta:
        contig2sequence = defaultdict(list)
        S_STRING = "S {} {}\n"
        L_STRING = "L {} + {} + 0M RC:i:{}\n"

        lengths = []

        with open(scaffolds_file) as scaffolds:
            header = scaffolds.readline().strip()[1:]
            #             headers.append(header)

            for line in tqdm(scaffolds, desc=f"Reading {scaffolds_file}"):
                if line.startswith(">"):
                    contig2sequence[header] = ''.join(contig2sequence[header])

                    if header in required_contigs_set:
                        save_gfa.write(S_STRING.format(*(header, contig2sequence[header])))
                        lengths.append(len(contig2sequence[header]))
                        save_fasta.write(f">{header}\n{contig2sequence[header]}\n")

                    header = line.strip()[1:]
                else:
                    contig2sequence[header].append(line.strip())
            else:
                contig2sequence[header] = ''.join(contig2sequence[header])
                if header in required_contigs_set:
                    lengths.append(len(contig2sequence[header]))
                    save_gfa.write(S_STRING.format(*(header, contig2sequence[header])))
                    save_fasta.write(f">{header}\n{contig2sequence[header]}\n")

        for _, row in tqdm(hic_data.data.iterrows(), total=hic_data.data.shape[0],
                           desc="Extracting Hi-C scores between contigs"):
            contig_1_2_score = row[[hic_data.node_1, hic_data.node_2, hic_data.score_column]]

            save_gfa.write(L_STRING.format(*contig_1_2_score))

    logger.info(f"Assembly contigs were saved to {save_fasta_file}")
    logger.info(f"Assembly graph was saved to {save_gfa_file}")

    if mimic:
        mimic_jgi(contigs_with_info_in_names=required_contigs,
                  # lengths=lengths,
                  save_name=os.path.join(saving_dir, "assembly_depth.tsv"),
                  initial_assembly_graph=initial_assembly_graph,
                  # contact_map=contact_map
                  )
    return contig2sequence


def abundance2jgi(abundances_file, save_dir, save_name="assembly_depth.txt"):
    save_file = os.path.join(save_dir, save_name)
    abundance = pd.read_csv(abundances_file, sep="\t", header=None, index_col=0)

    abundance.index.name = "contigName"
    abundance.columns = [f"file_{i}" for i in range(1, abundance.shape[1] + 1)]

    len_pattern = re.compile(r"length_(\d+)")
    regexp_vectorized = np.vectorize(lambda x: len_pattern.findall(x)[0])

    class RegexpError(Exception):
        def __init__(self):
            self.message = "Incorrect contig names as re couldn`t find length & coverage match"
            super().__init__(self.message)

    try:
        lengths = regexp_vectorized(abundance.index)

    except IndexError:  # re.findall has empty list <==> incorrect contig names
        raise RegexpError

    else:
        total_avg_depth = np.mean(abundance.values.astype(float), axis=1)
        abundance.insert(0, "contigLen", lengths)
        abundance.insert(1, "totalAvgDepth", total_avg_depth)
        abundance.to_csv(save_file, sep="\t", index=True)
        logger.info(f"Abundances were transformed to jgi format and saved to {save_file}")


def create_jgi_from_depth_file(depth_file, contigs_ordered, save_dir, save_name="assembly_depth.txt"):
    save_file = os.path.join(save_dir, save_name)
    depths = pd.read_csv(depth_file, delimiter="\t", index_col=0).iloc[:, :2]
    depths = depths.loc[contigs_ordered]
    depths.index.name = "contigName"
    depths.columns = ["contigLen", "totalAvgDepth"]
    depths["nodes.bam"] = depths["totalAvgDepth"]

    try:
        depths = depths.loc[contigs_ordered]
        depths.to_csv(save_file, sep="\t", index=True)
        logger.info(f"Depth file was read and jgi file was transformed with respect to it in {save_file}")

    except IndexError:
        breakpoint()


def plot_contigs_lengths_distribution(all_lengths, lengths_of_hic_scaffolds, outdir):
    def count_contigs_under_threshold(threshold, overall_space, limited_space):
        columns = ["Overall contigs", "Contigs without Hi-C links"]
        index = [" >= Threshold", " < Threshold"]

        over_ge = np.sum(overall_space >= threshold)
        over_l = np.sum(overall_space < threshold)

        lim_ge = np.sum(limited_space >= threshold)
        lim_l = np.sum(limited_space < threshold)

        data_report = pd.DataFrame([[over_ge, lim_ge], [over_l, lim_l]],
                                   columns=columns, index=index)

        return data_report

    count_t = partial(count_contigs_under_threshold,
                      overall_space=all_lengths,
                      limited_space=lengths_of_hic_scaffolds)

    thresholds = np.concatenate((np.arange(100, 600, 100), np.arange(500, 2500, 500), np.arange(3000, 11000, 1000),
                                 np.arange(20000, 5000000, 10000)))
    index_of_max = bisect.bisect_right(thresholds, np.max(all_lengths))
    thresholds = thresholds[:index_of_max]

    less = []

    logger.info(f"Estimating contigs distribution")
    for t in thresholds:
        less.append(list(count_t(t).iloc[1]))

    less = np.array(less)
    less = pd.DataFrame(columns=["Overall set", "With HiC links"], data=less, index=thresholds)
    less.index.name = "Threshold"

    less_then_3_sigma_border = np.sum(all_lengths < np.median(all_lengths) + 3 * np.std(all_lengths))

    fig = px.scatter(less, title="Distribution of # contigs with length less then threshold", log_x=True)
    fig.add_hline(y=less_then_3_sigma_border,
                  annotation_text=f"# contigs < median + 3 sigma of length ({less_then_3_sigma_border})",
                  annotation_position="top left")
    fig.update_layout(yaxis={"title": "Amount"}, legend={"title": "Contigs"}, font=dict(
        family="Proxima Nova",
        size=20,
        color="Dark Blue"
    ))
    fig.update_traces(mode='markers', marker_line_width=1, marker_size=15
                      )
    fig.write_image(os.path.join(outdir, "contig_lengths.jpg"))
