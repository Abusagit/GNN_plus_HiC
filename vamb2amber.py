#!/usr/bin/env python
# -*- coding: utf-8 -*-


import argparse
import pandas as pd
import os
import numpy as np
import re


class OutFileError(Exception):
    pass


class AmberProcessingError(Exception):
    pass


class AmberDataPreprocessor:
    def __init__(self, data, header_file):
        self.header = self.store_header(header_file)
        self.data = data

    @staticmethod
    def store_header(file):
        with open(file) as f_read:
            return ''.join(f_read.readline() for _ in range(2))

    def create_output_format(self, output_name):
        if output_name[-4:] != ".tsv":
            raise OutFileError("Incorrectly specified name (must be .tsv format)")

        self.data.to_csv(f"{output_name}", sep="\t", index=False)

        with open(f"{output_name}", "r+") as f:
            content = f.read()
            f.seek(0, 0)
            f.write(self.header.rstrip("\r\n") + "\n" + content)


def len_pattern_extractor(names: pd.Series):
    len_pattern = re.compile(r"length_(\d+)")

    regexp_vectorized = np.vectorize(lambda x: len_pattern.findall(x)[0])

    return regexp_vectorized(names)


def process_lengths_file(names: pd.Series):
    raise NotImplementedError


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str)
    parser.add_argument("-g", "--golden_standard", type=str)
    parser.add_argument("-o", "--outdir", type=str, default=os.path.join("./", "vamb_output.tsv"))
    parser.add_argument("-l", "--lengths", type=str, required=False)

    args = parser.parse_args()

    print("Starting transforming...")
    vamb_output = pd.read_csv(args.input, header=None, delimiter="\t")
    vamb_output[[0, 1]] = vamb_output[[1, 0]]
    vamb_output.rename(columns={0: "@@SEQUENCEID", 1: "BINID"}, inplace=True)

    if not args.lengths:
        print("No lengths given, trying to reach them from contig names...")

        try:
            lengths = len_pattern_extractor(vamb_output["@@SEQUENCEID"])

            exit_ = 1 / len(lengths)

        except ZeroDivisionError:
            raise AmberProcessingError(
                """
            Lengths weren`t found in contig names, tranfsormation is impossible with this set of data
            Prvode additional data for lengths of contigs
            """)
    else:
        lengths = process_lengths_file(vamb_output["@@SEQUENCEID"])

    vamb_output["_LENGTH"] = lengths

    print("Getting ready to the file writing...")
    amber_vamb_preproc = AmberDataPreprocessor(header_file=args.golden_standard, data=vamb_output)

    print(f"Writing at {args.outdir}")
    amber_vamb_preproc.create_output_format(args.outdir)


if __name__ == '__main__':
    main()
