#!/usr/bin/env python3
import argparse
import operator
import random

def parse_args():
    parser = argparse.ArgumentParser(description='Reduce 4D exons from BED file')
    parser.add_argument('-i', '--inputFile', required=True, help='Input BED file of 4D sites')
    parser.add_argument('-o', '--outputFile', help='Output file path')
    parser.add_argument('-n', '--numExons', type=int, default=100000, help='Number of exons to include')
    return parser.parse_args()

def reduce_4d_exons(input_file, output_file, num_exons):
    with open(input_file) as f:
        lines = f.readlines()

    random.shuffle(lines)
    keep_lines = sorted(lines[:num_exons], key=operator.itemgetter(0, 1))

    with open(output_file, 'w') as f:
        f.writelines(keep_lines)

def main():
    args = parse_args()
    output_file = args.outputFile or f'{args.inputFile.split(".")[0]}_reduced.bed'
    reduce_4d_exons(args.inputFile, output_file, args.numExons)

if __name__ == '__main__':
    main()
