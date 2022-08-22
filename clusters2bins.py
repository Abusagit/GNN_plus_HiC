import argparse
from collections import defaultdict
from tqdm import tqdm
# import os
from pathlib import Path


parser = argparse.ArgumentParser()

parser.add_argument("-c", "--clusters", help="File with binning results", type=str)
parser.add_argument("-f", "--fasta", help="File contatining contigs for bins (assembled contigs/scaffolds", type=str)

parser.add_argument("-o", "--out", type=str, help="Directory with bins: <out>/bins/", default=Path.cwd())

parser.add_argument("-l", "--mincontig", type=str, help="Minimal contig length value", default=2000)
parser.add_argument("-m", "--minbin", type=int, default=500_000)

args = parser.parse_args()


Path(args.out, "bins/").mkdir(parents=True, exist_ok=True)

fasta_nlines = sum(1 for _ in open(args.fasta, "rb"))
header2seq = {}

with open(args.fasta) as fasta_file:
    header = fasta_file.readline()

    tmp_sequence = []
    current_length = 0
    for line in tqdm(fasta_file, total=fasta_nlines-1, desc=f"Filtering contigs by length {args.mincontig}"):
        if line.startswith('>'):
            if current_length >= args.mincontig:
                header2seq[header.strip()[1:]] = ''.join([header] + tmp_sequence)
            
            tmp_sequence = []
            header = line
            current_length = 0

        else:
            tmp_sequence.append(line)
            current_length += len(line) - 1 # -1 for \n symbol
    else:
        if current_length >= args.mincontig:
            header2seq[header] = ''.join(tmp_sequence)



bin2seqs = defaultdict(list)
clusters_nlines = sum(1 for _ in open(args.clusters, "rb"))

with open(args.clusters) as clusters:
    for line in tqdm(clusters, total=clusters_nlines, desc="Managing clusters..."):
        seq, bin_, *_ = line.strip().split('\t')
        if seq in header2seq:
            bin2seqs[bin_].append(seq)


print(f"Filtering clusters by size {args.minbin}")

filtered_bins = filter(lambda bin_: sum(len(header2seq[c]) for c in bin2seqs[bin_]) >= args.minbin, bin2seqs.keys())

print(f"Writing filtered bins to {Path(args.out, 'bins/')}")

for i, filtered_bin in tqdm(enumerate(filtered_bins, start=1)):
    with open(Path(args.out, "bins/", f"{filtered_bin}.fa"), "w") as bin_file:
        bin_file.write(''.join(header2seq[contig] for contig in bin2seqs[filtered_bin]))

overall_bins = len(bin2seqs)
filtered_bins = i

print(f"Bins are saved at {Path(args.out, 'bins/')}")
print(f"{overall_bins=}, {filtered_bins=}")