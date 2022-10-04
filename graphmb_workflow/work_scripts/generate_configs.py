import argparse
from pathlib import Path

def get_parser():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-o", "--outdir", help="Directory for contigs")
    parser.add_argument("-i")
    
    return parser



def generate():
    pass