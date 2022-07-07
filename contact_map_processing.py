__doc__ = """Processing contact map and something for DMoN if needed"""

import numpy as np
import pandas as pd
import scipy.sparse as sparse
import os
import re

from tqdm import tqdm
from scipy.sparse import base
from collections import defaultdict
from sklearn.preprocessing import LabelEncoder
import subprocess

import logging

logger = logging.getLogger(__name__)


class ContactMap:
    DELIMITER = {"tsv": "\t", "csv": ","}
    SCALING_FUNC = {"log": np.log, "sqrt": np.sqrt}

    def __init__(self, path, feature_columns=None, node_1_column="FirstName",
                 node_2_column="SecondName",
                 score_column="SpadesScore",
                 scaling_method="log"):
        self.path = path

        self.node_1 = node_1_column
        self.node_2 = node_2_column

        self.feature_columns = feature_columns or []
        self.score_column = score_column
        self.scaling_func = scaling_method

        self.contigs_with_hic = None
        self.SHAPE = None
        self.nuinque_hic_contigs = None

        self.encoder = LabelEncoder()

        self.data = pd.read_csv(path, delimiter=ContactMap.DELIMITER[path.split(".")[-1]])

        logger.info(f"Scaling {self.score_column} by {self.scaling_func} scaling function")

        self.data[self.score_column] = self.data[self.score_column].apply(ContactMap.SCALING_FUNC[self.scaling_func])

        self.fit_encoder(node_1_column, node_2_column)

    def fit_encoder(self, node_1, node_2):
        logger.info("Fitting encoder...")
        unique_features = set(self.data[node_1].unique()) | set(self.data[node_2].unique())

        # fitting encoder:
        self.encoder.fit(list(unique_features))

        unique_features_amount = len(unique_features)
        self.nuinque_hic_contigs = unique_features_amount
        self.contigs_with_hic = unique_features
        self.SHAPE = (unique_features_amount, unique_features_amount)

        logger.info(f"Fitting complete! Found {self.nuinque_hic_contigs} contigs with Hi-C links")

    def get_sparse_adjacency_feature_matrices(self, threshold=0):
        """
        :param threshold: value SpadesScore must be greater or equal to for
        :return: sparse adjacency and sparse feature matrux for DMoN
        """

        h = len(set(self.data[self.node_1].unique()) | set(self.data[self.node_2].unique()))
        shape = (h, h)
        sparse_adjacency = sparse.lil_matrix(np.zeros(shape))  # !!!!!!!!!!!!!

        sparse_feature_matrix = sparse.lil_matrix(
            np.zeros((self.nuinque_hic_contigs, max(len(self.feature_columns) // 2, 1))))

        self.data[self.node_1] = self.encoder.transform(self.data[self.node_1].values)
        self.data[self.node_2] = self.encoder.transform(self.data[self.node_2].values)

        visited_pairs = set()
        for _, row in tqdm(self.data.iterrows(), total=self.data.shape[0]):

            node_1_id = row[self.node_1]
            node_2_id = row[self.node_2]

            if visited_pairs & {(node_1_id, node_2_id), (node_2_id, node_1_id)}:
                continue

            visited_pairs |= {(node_1_id, node_2_id), (node_2_id, node_1_id)}

            sparse_adjacency[node_1_id, node_2_id] = sparse_adjacency[node_2_id, node_1_id] = int(
                row[self.score_column] >= threshold)

            for feature_index, feature in enumerate(self.feature_columns[::2]):  # First node
                sparse_feature_matrix[node_1_id, feature_index] = row[feature]

            for feature_index, feature in enumerate(self.feature_columns[1::2]):  # Second node
                sparse_feature_matrix[node_2_id, feature_index] = row[feature]

        sparse_adjacency = sparse.csr_matrix(sparse_adjacency)
        sparse_feature_matrix = sparse.csr_matrix(sparse_feature_matrix)

        return sparse_adjacency, sparse_feature_matrix
