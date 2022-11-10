from collections import defaultdict
import pandas as pd
from pathlib import Path
import argparse
from itertools import takewhile, repeat
from tqdm import tqdm

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta")
    parser.add_argument("-c", "--clusters")
    parser.add_argument("--outdir", default=Path(Path.cwd(), "bins"), type=Path)
    parser.add_argument("--minbin", default=500_000, type=int)
    
    return parser


def linecount(filename):
    print(f"Counting lines from {filename}")
    with open(filename, "rb") as file:
        bufgen = takewhile(lambda x: x, (file.raw.read(1024*1024) for _ in repeat(None)))
    
        return sum( buf.count(b'\n') for buf in bufgen )

def read_fasta(file):
    sequences = dict()
    
    
    lines = linecount(file)
    
    with open(file) as fasta_read:
        header = fasta_read.readline().strip()[1:]
        _seq = []
        
        lengths = dict()
        for line in tqdm(fasta_read, total=lines-1, desc=f"reading {file}..."):
            if line[0] == '>':
                sequence = ''.join(_seq)

                sequences[header] = sequence
                lengths[header] = len(sequence)
                sequence = ''
                _seq = []
                
                header = line.strip()[1:]
            else:
                _seq.append(line.strip())
        else:
            sequence = ''.join(_seq)
            sequences[header] = sequence
            lengths[header] = len(sequence)
            
        del _seq, header
    
    return sequences, lengths

# Path().mkdir()
def write_bins(fasta_dict, clusters_dict, contiglen, minbin_len, outdir):
    outdir.mkdir(parents=True, exist_ok=False)
    
    
    for bin_, contigs in tqdm(clusters_dict.items(), total=len(clusters_dict), desc=f"Writing bins to {str(outdir)}"):
        if sum(contiglen[c] for c in contigs) < minbin_len:
            continue
        
#        breakpoint()
        with open(outdir / f"{bin_}.fa", "w") as f_write:
            f_write.write('\n'.join(">{}\n{}".format(contig, fasta_dict[contig]) for contig in contigs))
    
    


def main():
    
    parser = get_parser()
    args = parser.parse_args()
    
    clusters = pd.read_csv(args.clusters, index_col=0, header=None, sep='\t')
    contig2bin = dict(zip(clusters.index, clusters.iloc[:, 0]))
    
    bin2contig = defaultdict(list)
    for contig, bin_ in contig2bin.items():
        bin2contig[bin_].append(contig)
    
    del clusters
    
    fasta, contiglen = read_fasta(args.fasta)
    
    write_bins(fasta, bin2contig, contiglen, args.minbin, args.outdir)
    
    
if __name__ == "__main__":
    main()