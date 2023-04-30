#!/usr/bin/env python3

import sys
import os
import time
import random
import numpy
import gzip
import math
import pandas
import argparse

"""
This program converts a .gff or .gff.gz file into a .bed file.
"""

def parse_args():
    """
    Parses command-line arguments.

    Returns:
        args: A Namespace object containing the command-line arguments.
    """
    parser = argparse.ArgumentParser(description='Convert GTF annotation file to BED6 format for use by HAL to extract 4d sites.')
    parser.add_argument('-i', '--input_file', required=True, help='Input GTF file, compressed or otherwise.')
    parser.add_argument('-o', '--output_file', help='Output file path to .bed file created by this function.')
    parser.add_argument('-l', '--loci', type=int, default=100000, help='Number of loci to include for training phyloP model.')
    return parser.parse_args()

def convert_file(input_file, output_file, loci):
    """
    Converts a .gff or .gff.gz file to a .bed file.

    Args:
        input_file: A string containing the path to the input file.
        output_file: A string containing the path to the output file.
        loci: An integer specifying the number of loci to include for training phyloP model.
    """

    # We only require a subset of the 4d sites
    # For simplicity, get all CDSs in frame 0:
    # Names taken from https://useast.ensembl.org/info/website/upload/gff.html
    df = pandas.read_csv(input_file,sep='\t',names=['seqname','source','feature','start','end','score',
        'strand','frame','attribute'])
    df = df.loc[(df['feature'] == 'CDS') & (df['frame'] == '0')]
    # GTFs are 1-based, and .beds are 0-based: convert:
    df['start'] = df['start']-1
    df['end'] = df['end']-1
    df['blockSizes'] = df['end']-df['start']
    df = df.reset_index(drop=True)
    df['uid'] = ['CDS_'+str(x) for x in df.index]
    df['thickStart'] = df['start']
    df['thickEnd'] = df['end']
    df['itemRgb'] = '255,0,0'
    df['blockCount'] = 1
    df['blockStarts'] = 0
    df = df.sort_values(by=['blockSizes'],ascending=False).head(loci)
    df[['seqname','start','end','uid','score','strand','blockCount','blockSizes','blockStarts']].to_csv(output_file,sep='\t',index=False,header=None)

def main():
    args = parse_args()
    convert_file(args.input_file, args.output_file, args.loci)

if __name__ == '__main__':
    main()