#!/usr/bin/env python3
import argparse


class CommandLine:
    """Handles the input arguments from the command line. Manages 
    the argument parser.
    Methods:
    Other than initialization, no methods are present, as its purpose is 
    simply to handle what is passed into the command line and pass that 
    into the class that performs the searching algorithm."""

    def __init__(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("-s", "--ssFile", help="Input .ss of high-confidence tRNAs from tRNAscan-SE 2.0.")
        parser.add_argument("-b", "--bedFile", help="Input .bed of high-confidence tRNAs from tRNAscan-SE 2.0.")
        parser.add_argument("-o", "--outFile", help="File to be analyzed by RNAfold to get MFE data for tRNAs.", default='tRNAFoldIn.txt')
        self.args = parser.parse_args()


class FileConverter:
    """
    Primary class where filetype is converted.
    """

    def __init__(self, ssFile, bedFile, outFile):
        self.ssFile = ssFile
        self.bedFile = bedFile
        self.outFile = outFile

    def process_str(self, my_string):
        return my_string.replace('>', '(').replace('.', '.').replace('<', ')')

    def convert_file(self):
        coord_to_tRNA = {}
        with open(self.bedFile) as f:
            for line in f:
                split_line = line.strip().split('\t')
                for k in range(int(split_line[1]), int(split_line[2]) + 1):
                    coord_to_tRNA[f'{split_line[0]}__{k}'] = str(split_line[3])
        out_string = ''
        with open(self.ssFile) as f:
            tRNA = ''
            for line in f:
                split_line = line.strip().split()
                if len(split_line) > 2 and 'Length' in split_line[2]:
                    chrom = split_line[0].split('.trna')[0]
                    coords = ((split_line[1])[1:-1]).split('-')
                    start, end = min(int(coords[0]), int(coords[1]))
                    if f'{chrom}__{end + 25}' in coord_to_tRNA:
                        tRNA = coord_to_tRNA[f'{chrom}__{end + 25}']
                if len(tRNA) > 1 and split_line[0].startswith('Seq:'):
                    seq = split_line[1].upper()
                if len(tRNA) > 1 and split_line[0].startswith('Str:'):
                    out_string += f'>{tRNA}\n{seq}\n{self.process_str(split_line[1])}\n'
                    tRNA = ''
        with open(self.outFile, 'w') as f:
            f.write(out_string)


def main():
    """
    Initializes a CommandLine object and passes the provided 
    arguments into a new FileConverter object and calls main method.
    """
    command_line = CommandLine()
    file_converter = FileConverter(command_line.args.ssFile, command_line.args.bedFile, command_line.args.outFile)
    file_converter.convert_file()


if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main()
