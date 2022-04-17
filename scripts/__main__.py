
from functions import *
import argparse
# What is this files? https://stackoverflow.com/questions/4042905/what-is-main-py

################################################################################
## Command line arguments
################################################################################
def is_file(string):
    """Check if the input file provided is a file
    """
    import os

    # Check the extension of the file
    if string.endswith("fa") or string.endswith("fasta"):
        return string
    else:
        msg = f"{string} is has not a valid input_file format (.fa/.fasta)"
        raise arg.parse.ArgumentTypeError(msg)

    if not os.path.isfile(string):
        msg = f"{string} is not an existing input file"
        raise arg.parse.ArgumentTypeError(msg)

parser = argparse.ArgumentParser(description=
    "Flexign provides a flexibility score given a protein sequence or a uniprot identifier. It can also take a protein family from which the consensus sequence is retrieved and used as input")

# Input argument
parser.add_argument('-i', '--input', required = True,
                    type = is_file,
                    dest = "infile",
                    action = "store",
                    default = None,
                    help = "Input a FASTA file")

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

input_file = options.infile
output_file = options.outfile

################################################################################
## MAIN CODE
################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get bfactors for target protein
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Get Homologs
protein_codes, protein_chains, hits_dict = get_pdb_homologs(input_file)
# Download homologs
pdb_outfiles = download_pdb_files(protein_codes)
# Apply resolution quality filter
final_protein_codes, final_chains = pdb_quality_filter(pdb_outfiles, protein_codes, hits_dict)
# Download filtered homologs
pdb_outfiles = [f"structures/pdb{code}.ent" for code in final_protein_codes]
# Get fasta sequences for multiple sesquence alignment
pdb_seqs_filename = get_pdb_sequences(final_protein_codes, final_chains)
# Perfom Multiple Sequence Aligment
msa_obj= msa(pdb_seqs_filename)
# Obtain position of identifical residues & the length of the pdb residues
all_sim, all_con, all_rest, all_len_pdb, all_seqs = clustal_annotations(msa_obj)

# Obtain the mean of normalized b-factors for residues in regions with similarities
# Get the b-factors for all residues and normalize

all_bfactors = [get_bfactors(prot[0], prot[1], prot[2], prot[3]) for prot in
                 list(zip(pdb_outfiles,final_protein_codes, final_chains, all_len_pdb))]
all_norm_bfactors = [normalize_bfactor(bfactors) for bfactors in all_bfactors]

## Get the b-factor only for regions of similarity
final_bfactors = get_modified_bfactors(all_norm_bfactors, all_sim, all_con, all_rest, all_seqs)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get the flex score
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## GEt get_residues
residues = [get_bfactors2(prot[0], prot[1], prot[2], prot[3], prot[4]) for prot in
                 list(zip(pdb_outfiles,final_protein_codes, final_chains, all_sim, all_len_pdb))]

## Compute F score validation
r0 = []
for i in residues[0]:
    r0.append(i.get_resname())

F5_val = flexibility_val(r0, window_size=5)
F1_cval = flexibility_val(r0, window_size=1)

## Compute the F score
F5,flex_scores = flexibility(final_bfactors, window_size = 5)
F3 = flexibility(final_bfactors, window_size = 3)
F1 = flexibility(final_bfactors, window_size = 1)

## F scores normalized
norm_flex_scores = scale_function(flex_scores)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get Secondary Structure
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from pathlib import Path
from glob import glob
# Change the extension of pdb '.ent' to '.pdb'
files = glob("structures/*.ent")
for name in files:
    p = Path(name)
    s = p.rename(p.with_suffix('.pdb'))

pdb_outfiles = [f"structures/pdb{code}.pdb" for code in final_protein_codes]


sstructures = [get_sstructure(prot[0], prot[1], prot[2]) for prot in list(zip(pdb_outfiles, final_protein_codes, final_chains))]
# Get Secondary structure of regions of similarity
all_sim_sstructures = list()
for val in list(zip(sstructures, all_sim)):
    sstructure = val[0]
    similarities = val[1]

    sim_sstructure = [sstructure[i] for i in similarities]
    all_sim_sstructures.append(sim_sstructure)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get Hidrobicity
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hydroph_scores_aa, gravy = get_hydrophobicity(input_file)
