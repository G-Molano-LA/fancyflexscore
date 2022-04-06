#!/usr/bin/env python3

## Search candidates that has similar structure with our target protein

# Python modules
import logging
import argparse
import os
from structural_py.lib.Bio import SeqIO

############################################################
## Command line arguments
############################################################

parser = argparse.ArgumentParser(description=
    "Flexign provides a flexibility score given a protein sequence or a uniprot identifier. It can also take a protein family from which the consensus sequence is retrieved and used as input")

# Input argument
parser.add_argument('-i', '--input', required = True,
                    dest = "infile",
                    action = "store",
                    default = None, 
                    help = "Input the FASTA file")

# Output argument
parser.add_argument('-o', '--output', required = True,
                    dest = "outfile",
                    action = "store",
                    default = "flexign_score_output.txt", 
                    help = "Output file with the flexibility score and a graphic representation of the flexibility by aminoacid")

# Verbose argument
parser.add_argument('-v', '--verbose',
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Select this option to have a follow up of the program")

options = parser.parse_args()

############################################################
## Accessing the files
############################################################

# Creating input and output variables
input_file = options.infile
output_file = options.outfile

# Checking if the input is a file
if input_file:
    if os.path.isfile(input_file):
        fasta_sequence = SeqIO.parse(open(input_file),'fasta')

