import pandas as pd
import argparse
import numpy as np



def get_parser():
    parser = argparse.ArgumentParser("raw bin3C contact map -> uniform format")
    
    parser.add_argument("-i", "--input_contact_map")
    parser.add_argument("-d", "--depths")
    parser.add_argument("-o", "--out", default="bin3c_contact_map.tsv")
    
    return parser

def main():
    parser = get_parser()
    args = parser.parse_args()
    
    bin3c_contact_map = pd.read_csv(args.input_contact_map, header=None, sep='\t')
    contig_depths_lengths = pd.read_csv(args.depths, index_col=0, sep='\t')
    
    # breakpoint()
    
    contigs_1 = contig_depths_lengths.loc[bin3c_contact_map.iloc[:, 0].values, ["contigLen", "totalAvgDepth"]]
    contigs_2 = contig_depths_lengths.loc[bin3c_contact_map.iloc[:, 1].values, ["contigLen", "totalAvgDepth"]]
    
    length_1, cov_1 = contigs_1["contigLen"].values.astype(int), contigs_1["totalAvgDepth"].values
    length_2, cov_2 = contigs_2["contigLen"].values.astype(int), contigs_2["totalAvgDepth"].values
    scores = bin3c_contact_map.iloc[:, 2]
    new_contact_map = pd.DataFrame(data={"FirstName": bin3c_contact_map.iloc[:, 0].values,
                                         "SecondName": bin3c_contact_map.iloc[:, 1].values,
                                         "FirstCoverage": cov_1,
                                         "SecondCoverage": cov_2,
                                         "FirstLength": length_1,
                                         "SecondLength": length_2,
                                         "SpadesScore": scores})
    
    print(f"Saving at {args.out}...")
    new_contact_map.to_csv(args.out, index=False, header=True, sep='\t')
    return 0

if __name__ == "__main__":
    main()