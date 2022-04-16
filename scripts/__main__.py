from functions import *

# What is this files? https://stackoverflow.com/questions/4042905/what-is-main-py


protein_codes = ["1xb7", "2ewp", "3d24"] # PDB codes in lowercasese
protein_chains = ["A", "E", "A"]
# Output format of pdb_files
pdb_outfiles = [f"structures/pdb{code}.ent" for code in protein_codes]

# Get PDB files
pdb_seqs_filename = get_pdb_sequences(protein_codes, protein_chains, pdb_outfiles)
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
