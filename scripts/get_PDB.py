# Obtain the PDB model of the target protein
## Obtain PDB from alpha-fold is the sequence is not known
## Obtain PDB from database if the protein is known

# First of all
# pip install biopython

# Modules
## To track events
import logging
## PDB files
from Bio.PDB.PDBParser import PDBParser         # read PDB files
from Bio.PDB.PDBList import download_pdb_files   # to get PDB file from PDB server
## To download files from web
import requests

# General configuration of logging module. More info: https://docs.python.org/3/howto/logging.html#logging-basic-tutorial
logging.basicConfig(filename='flexscore.log', encoding = 'utf-8', level=logging.info)

################################################################################
## Here starts code
################################################################################

# Inputs
candidates = ["candidate1", "candidate2"]
target     = ["target"]



# 1. Obtain PDB files for candidates and target
logging.info("Obtaining PDB files for candidates and target proteins")

# Outputs
## ????
pdb_files = [ name for name in candidates] # mejor hacer una copia o algo asi?
pdb_files.insert(1, target)

# 1.1. Obtain PDB file for candidate proteins
logging.info("Obtaining PDB files for candidates")
## Save PDB results in a dictionary
pdb_structures = {}
## Get PDBs
r = PDBList() # Object instance
r.download_pdb_files(pdb_codes = pdb_files, file_format = "pdb")
## ?? Check here what happens when no PDB
## ?? check the directory created. If not like, use pdir argument


# 1.2. Obtain PDB file for target protein
logging.info("Obtaining PDB files for target protein")

## If there is not pdb file, get from alpha fold
if target_pdb is empty:
    logging.info("Target protein has not PDB file. Getting structure from AlphaFold database")

    # Obtain pdb file from AlphaFold
    ## Define the remote file to retrieve
    url = f"https://alphafold.ebi.ac.uk/files/AF-{target}-F1-model_v2.pdb"
    # Define the local filename to save data
    local_file = f'{target}.pdb'
    # Make http request for remote file data
    target_data = requests.get(url)
    # Save file data to local copy
    with open(local_file, 'wb')as file:
        file.write(target_data.content)



# 2. Read PDB files
logging.info("Reading PDB files")
## PBD parser instance
p = PDBParser(PERMISSIVE = TRUE)


# Parser PDB files
for i in range(0, len(candidates), 1):
    pdb_structures[candidates[i]] = p.get_structure(candidates[i], pdb_files[i])
