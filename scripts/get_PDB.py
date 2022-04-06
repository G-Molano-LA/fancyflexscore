# Obtain the PDB model of the target protein
## Obtain PDB from alpha-fold if the sequence is not known
## Obtain PDB from database if the protein is known

# First of all
# pip install biopython

# Modules
## To track events
import logging
## PDB files
from Bio.PDB.PDBParser import PDBParser         # read PDB files
from Bio.PDB.PDBList import PDBList   # to get PDB file from PDB server
## Uniprot codes
import Bio.SwissProt as sp
## To download files from web
import requests
## Catch the stdout
from io import StringIO
import sys
## Regex
import re

# class capturing_output(list):
#     """Captures STDOUT by creating a new STDOUT
#
#     The special methods (__enter__, __exit__) allows create class instances by
#     using the 'with' statement.
#
#     Example:
#         with capturing_output() as out:
#             do_something(my_object)
#         out.getvalue()
#     """
#
#     def __enter__(self):
#         # Save the current stdout in a class atttribtue
#         self._stdout = sys.stdout
#         # and create a new one, save it as a class attribute and replace the previous one
#         sys.stdout = self._stringio = StringIO()
#         return self
#     def __exit__(self, *args):
#         # return a list object
#         self.extend(self._stringio.getvalue().splitlines())
#         # Delete created STDOUT to free up some memory
#         del self._stringio
#         # Restore the initial stdout
#         sys.stdout = self._stdout


def get_stream_data(url, outfile_path):
    """Download data from websites
    """
    import requests

    # Make http request for remote file data
    target_data = requests.get(url)
    # Save file data to local copy
    with open(outfile_path, 'wb')as file:
        file.write(target_data.content)

# General configuration of logging module. More info: https://docs.python.org/3/howto/logging.html#logging-basic-tutorial
logging.basicConfig(filename='flexscore.log', encoding = 'utf-8', level=20) # level logging.info

################################################################################
## Here starts code
################################################################################

# Inputs
protein_codes = ["1A3N", "2DHB"] # PDB codes

# Obtain PDB code for target protein, if any
## Obtain uniprot file
# url       = f"https://www.uniprot.org/uniprot/{target}.txt"
# localfile = f"Uniprot/{target}.txt"
# get_stream_data(url, localfile)
# ## Read uniprot entrance
# with open(localfile) as fh:
#     uniprot_record = sp.read(fh)
#     print(uniprot_record.references)
#     print(uniprot_record.cross_references)


# 1. Obtain PDB files for candidates and target
logging.info("Obtaining PDB files for candidate and target proteins")

# Outputs
## ????
# pdb_files = [ name for name in candidates] # mejor hacer una copia o algo asi?
# pdb_files.insert(1, target)

# 1.1. Obtain PDB file for target and candidates proteins
logging.info("Obtaining PDB files for candidates")
## Save PDB results in a dictionary
pdb_structures = {}
## Get PDBs
r = PDBList() # Object instance
r.download_pdb_files(pdb_codes = protein_codes , file_format = "pdb",
                    pdir = "structures/", overwrite=True) # creates the directory if do not exists


## If there is not pdb file, get from alpha fold
# if "doesn't exists" in captured_out[1]:
#     """
#     The download_pdb_files function always tries to download the pdb file,
#     printing:
#         "Downloading PDB structure 'pdb_code'...""
#     We the structure is not found, prints:
#         "Desired structure doesn't exists"
#     Therefore, if there is no structure for the target, we will download the structure
#     from alphafold
#     """
#     # Get the protein_code
#     logging.info(f"{target} has not PDB file. Getting structure from AlphaFold database...")
#     # Obtain pdb file from AlphaFold
#     ## Define the remote file to retrieve
#     url = f"https://alphafold.ebi.ac.uk/files/AF-{target}-F1-model_v2.pdb"
#     # Define the local filename to save data
#     local_file = f'{target}.pdb'
#     # Make http request for remote file data
#     get_stream_data(url, localfile)



# 2. Read PDB files
logging.info("Reading PDB files")
## PBD parser instance
p = PDBParser(PERMISSIVE = TRUE)


# Parser PDB files
for i in range(0, len(candidates), 1):
    pdb_structures[candidates[i]] = p.get_structure(candidates[i]+, pdb_files[i])
