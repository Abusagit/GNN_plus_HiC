import argparse
import pickle
from pathlib import Path

from tqdm import tqdm

import scipy.sparse as sparse


# package_dir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
# sys.path.insert(0, os.path.join(package_dir, "working_dir"))
def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser("Converting hicbin encoded contact map to uniform format")
    
    parser.add_argument("-n", "--names")
    parser.add_argument("-m", "--contact-map", "--contact_map")
    parser.add_argument("-o", "--outname")
    
    return parser


def read_seqnames(seq_info_path) -> list:
    with open(seq_info_path, "rb") as handler:
        
        seq_names_length_coverages = pickle.load(handler)
        
    return seq_names_length_coverages
    
def read_contact_map(old_contact_map_path):
    
    matrix = sparse.load_npz(old_contact_map_path)
    rows, cols = matrix.nonzero()
    print(matrix.data)

    
    return rows, cols, matrix.tocsr()

def fill_contact_map(rows, cols, matrix, info, outname):
    name_splitted = outname.split('/')
    if len(name_splitted) > 1:
        output = Path(*name_splitted[:-1])
        output.mkdir(exist_ok=True, parents=True)
    else:
        output = Path.cwd()
        
    output /= name_splitted[-1]
    
    with open(output, "w") as contact_map_file:
        contact_map_file.write("FirstName\tSecondName\tFirstCoverage\tSecondCoverage\tFirstLength\tSecondLength\tSpadesScore\n")
        
        LINE_STR = '\t'.join(["{}", "{}", "{}", "{}", "{}", "{}", "{}\n"])
        for row, col in tqdm(zip(rows, cols), total=len(rows)):
            score = matrix[row, col]

            name_1, len_1, cov_1 = info[row]
            name_2, len_2, cov_2 = info[col]
            
            line = LINE_STR.format(name_1, name_2, cov_1, cov_2, len_1, len_2, score)
            breakpoint()
            
            contact_map_file.write(line)
    
            
            
    
    
def main():
    parser = get_parser()
    args = parser.parse_args()
    
    
    seqs_infos = read_seqnames(args.names)
    rows, cols, matrix = read_contact_map(args.contact_map)
    fill_contact_map(rows=rows, cols=cols, matrix=matrix,
                     info=seqs_infos, outname=args.outname)


main()