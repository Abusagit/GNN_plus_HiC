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

logger = logging.getLogger(__name__)


class AmberDataPreprocessor:
    def __init__(self, data, header_file):
        self.header = self.store_header(header_file)
        self.data = data

    @staticmethod
    def store_header(file):
        with open(file) as f_read:
            return ''.join(f_read.readline() for _ in range(2))

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
    print(f"File saved with name '{out_file_name}.npz'")


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


def mimic_jgi_get_depth_from_names(contigs_with_info_in_names, lengths=None, save_name="features.txt",
                                   write_option=True):
    """
    Only suitable for data with the specidic pattern
    """

    len_cov_pattern = re.compile(r"length_(\d+)_cov_(\d+\.\d*)")

    regexp_vectorized = np.vectorize(lambda x: len_cov_pattern.findall(x)[0])

    class RegexpError(Exception):
        def __init__(self):
            self.message = "Incorrect contig names as re couldn`t find length & coverage match"
            super().__init__(self.message)

    try:
        if not lengths:
            lengths, coverages = regexp_vectorized(contigs_with_info_in_names)
        else:
            _, coverages = regexp_vectorized(contigs_with_info_in_names)
    except IndexError:  # re.findall has empty list <==> incorrect contig names
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


def create_gfa_file(hic_data, required_contigs, scaffolds_file, saving_dir="./", mimic=True,
                    save_gfa_file="contact_graph.gfa",
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

            for line in tqdm(scaffolds):
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

        for _, row in tqdm(hic_data.data.iterrows(), total=hic_data.data.shape[0]):
            contig_1_2_score = row[[hic_data.node_1, hic_data.node_2, hic_data.score_column]]

            save_gfa.write(L_STRING.format(*contig_1_2_score))

    logger.info(f"Assembly contigs were saved to {save_fasta_file}")
    logger.info(f"Assembly graph was saved to {save_gfa_file}")

    if mimic:
        mimic_jgi_get_depth_from_names(contigs_with_info_in_names=required_contigs,
                                       lengths=lengths,
                                       save_name=os.path.join(saving_dir, "assembly_depth.tsv"),
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
    less = pd.DataFrame(columns=["Whole contigs set", "Without HiC links"], data=less, index=thresholds)
    less.index.name = "Threshold"

    less_then_3_sigma_border = np.sum(all_lengths < np.median(all_lengths) + 3 * np.std(all_lengths))

    fig = px.scatter(less, title="Distribution of # contigs with length less then threshold", log_x=True)
    fig.add_hline(y=less_then_3_sigma_border,
                  annotation_text=f"# contigs < median + 3 sigma of length ({less_then_3_sigma_border})",
                  annotation_position="top left")
    fig.update_layout(yaxis={"title": "Amount"}, legend={"title": "Contigs"})
    fig.write_image(os.path.join(outdir, "contig_lengths.jpg"))
