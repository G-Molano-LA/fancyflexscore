
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
# Get Homologs
protein_codes, protein_chains, hits_dict = get_pdb_homologs(input_file)
# Download homologs
pdb_outfiles = download_pdb_files(protein_codes)
# Apply resolution quality filter
final_protein_codes, final_chains = pdb_quality_filter(pdb_outfiles, protein_codes, hits_dict)
# Download filtered homologs
pdb_outfiles = download_pdb_files(final_protein_codes)

# Get fasta sequences for multiple sesquence alignment
pdb_seqs_filename = get_pdb_sequences(protein_codes, protein_chains)
# Perfom Multiple Sequence Aligment
msa_obj= msa(pdb_seqs_filename)
# Obtain position of identifical residues & the length of the pdb residues
all_sim, all_aln_res = clustal_annotations(msa_obj)

# Obtain the mean of normalized b-factors for residues in regions with similarities
# Get the b-factors for all residues and normalize

all_bfactors = [get_bfactors(prot[0], prot[1], prot[2], prot[3]) for prot in
                 list(zip(pdb_outfiles,protein_codes, protein_chains, all_aln_res))]
all_norm_bfactors = [normalize_bfactor(bfactors) for bfactors in all_bfactors]

## Get the b-factor only for regions of similarity
all_sim_bfactors = list()

for val in list(zip(all_norm_bfactors, all_sim)):
    norm_bfactors = val[0]
    similarities  = val[1]

    sim_bfactors = [norm_bfactors[i] for i in similarities]
    all_sim_bfactors.append(sim_bfactors)

# Compute the mean
final_bfactors = [sum(sim_bfactors) for sim_bfactors in zip(*all_sim_bfactors)]

# Get Secondary Structure
pdb_outfiles = [f"structures/pdb{code}.pdb" for code in protein_codes]
sstructures = [get_sstructure(prot[0], prot[1], prot[2]) for prot in list(zip(pdb_outfiles, protein_codes, protein_chains))]
# Get Secondary structure of regions of similarity
all_sim_sstructures = list()
for val in list(zip(sstructures, all_sim)):
    sstructure = val[0]
    similarities = val[1]

    sim_sstructure = [sstructure[i] for i in similarities]
    all_sim_sstructures.append(sim_sstructure)

# Get thee flex score
## GEt get_residues
residues = [get_bfactors2(prot[0], prot[1], prot[2], prot[3], prot[4]) for prot in
                 list(zip(pdb_outfiles,protein_codes, protein_chains, all_sim, all_aln_res))]

## Compute F score validation
r0 = []
for i in residues[0]:
    r0.append(i.get_resname())

F5_val = flexibility_val(r0, window_size=5)
F1_cval = flexibility_val(r0, window_size=1)

## Get the b-factor only for regions of similarity
all_sim_bfactors = list()

for val in list(zip(all_norm_bfactors, all_sim)):
    norm_bfactors = val[0]
    similarities  = val[1]

    sim_bfactors = [norm_bfactors[i] for i in similarities]
    all_sim_bfactors.append(sim_bfactors)

## Compute the F score
F5,flex_scores = flexibility(all_sim_bfactors, window_size = 5)
F3 = flexibility(all_sim_bfactors, window_size = 3)
F1 = flexibility(all_sim_bfactors, window_size = 1)

## F scores normalized
norm_flex_scores = scale_function(flex_scores)

# Get Hidrobicity
hydroph_scores_aa, gravy = get_hydrophobicity(input_file)
