#!/usr/bin/env python3

# Obtain regions of similarity:
## Structural MSA implementation

# Python modules
## To track events
import logging
## PDB biopyhton  (upload only the required ones)
import Bio.PDB
##
import argparse
import os
#from structural_py.lib.Bio import SeqIO
## Catch the stdout
from io import StringIO
import sys
## Regex
# mport re
# MSA
from caretta import multiple_alignment

multiple_alignment.StructureMultiple.align_from_pdb_files("1xb7.pdb")
############################################################
## Accessing the files
############################################################

# Creating input and output variables
# input_file = options.infile
# output_file = options.outfile

# Checking if the input is a file
# if input_file:
#     if os.path.isfile(input_file):

# PDB code
pdb_code = ["1xb7","2pjl","3d24"]
pdb_filename = ("%s.pdb,%s.pdb,%s.pdb" %(pdb_code[0],pdb_code[1],pdb_code[2]))
pdb_filename = pdb_filename.split(',')
print(pdb_filename)



# def superimposition(pdb_code,pdb_filename):
#
#     # Load PDB file
#     structure = Bio.PDB.PDBParser().get_structure(pdb_code, pdb_filename)
#     dic_of_pairs = {}
#     print(len(structure))
#     for model in structure:
#         input_pair_list = []
#         for chain in model:
#             input_pair_list.append(chain)
#         dic_of_pairs[model]= input_pair_list
#     return dic_of_pairs
#
# print(superimposition(pdb_code,pdb_filename))


######## NOTAS ########
# There are some methods (scripts like modeller) to perform the structural alignment:
# matt, caretta and tm-aling
# pdb files in lowecases
#
