
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
# Apply resolution quality filter
final_protein_codes, final_chains = pdb_quality_filter(protein_codes, hits_dict)
# Download filtered homologs
pdb_outfiles = [f"structures/pdb{code}.ent" for code in final_protein_codes]
# Get fasta sequences for multiple sesquence alignment
pdb_seqs_filename = get_aln_sequences(input_file, pdb_outfiles, final_protein_codes, final_chains)
# Perfom Multiple Sequence Aligment
msa_obj= msa(pdb_seqs_filename)
# Obtain sequences and len of pdb alignment
all_seqs, all_len_pdb, all_matrices = get_pdb_seq(msa_obj)
target_seq_len = all_len_pdb.pop(0)
# Get the b-factors for all residues and normalize
all_bfactors = [get_bfactors(prot[0], prot[1], prot[2], prot[3]) for prot in
                 list(zip(pdb_outfiles,final_protein_codes, final_chains, all_len_pdb))]
all_norm_bfactors = [normalize_bfactor(bfactors) for bfactors in all_bfactors]

# Get the final b-factors
msa_seqs = [obj.seq for obj in msa_obj]


target_norm_bfactors = [None]*(target_seq_len+1)

# Modify the bfactors according to msa
final_bfactors = get_modified_bfactors(msa_seqs, all_norm_bfactors, target_norm_bfactors, all_matrices)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get the flex score
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Compute the F score
F7,flex_scores = flexibility(final_bfactors, window_size =7)
F3 = flexibility(final_bfactors, window_size = 3)
F1 = flexibility(final_bfactors, window_size = 1)

## F scores normalized
norm_flex_scores = scale_function(flex_scores)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get Secondary Structure
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get pdb file from alpha fold
# !!! fasta id must contain only the uniprot code. If there is chain, mustbe specified by uniprotid_chain
target_id, target_seq = extract_fasta_sequence(input_file)
# Checking if the target file has chain or not
split_id = target_id.split('_')

if len(split_id) == 2:
    target_id = split_id[0]
    target_chain = split_id[1]
    # Get pdb from AlphaFold
    local_file = get_pdb_from_alphafold(target_id)
    # get sstructure
    sstructures = get_sstructure(local_file, target_id, target_chain)
else:
    # Get pdb from AlphaFold
    local_file = get_pdb_from_alphafold(target_id)
    # get sstructure
    sstructures = get_sstructure(local_file, target_id)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get Hidrobicity
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hydroph_scores_aa, gravy = get_hydrophobicity(input_file)
