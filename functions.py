def get_pdb_homologs(input_file, E_VALUE_THRESH = 1e-20):
    """Obtain PDB codes and chains from protein homologs of the query protein

    This function conducts a BLASTP against the PDB database locally and retrieves the unique
    hits with their respective chains. Proteins with homology in several chains are not considered unique
    and only the first chain is kept.
    In case less than 7 hits, with a lower threshold than 1e-20, are found two psi-blasts are conducted. First,
    against the Uniprot database to obtain the 3rd iteration PSSM and then against the PDB database using the
    obtained PSSM. This is only conducted to fill the list with remote homologs in case not enough homologs are obtained
    with the BLASTP.

    Input: FASTA file with the query protein sequence

    Return: List of protein codes, list of protein chains of the homologs and dictionary with the proteins as keys and
    chains as values
    """

    # Python modules
    import os
    import glob # for matching wildcards
    import shutil
    import subprocess
    from Bio import SeqIO
    from Bio.Blast.Applications import NcbiblastpCommandline
    from Bio.Blast import NCBIXML
    from Bio.Blast.Applications import NcbipsiblastCommandline

    # Checking if the input is a file or a sequence:
    if os.path.isfile(input_file):
        # Get BLASTP db
        if os.path.isdir('./db/blastp'):
            if os.listdir("./db/blastp"):
                print("Using the existing blastp database located at './db/blastp'")
            else:
                raise OSError("Not database found! Please delete the ./db/blastp folder to allow the creation of a new one")
        else:
            subprocess.run(["wget", "https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz"], check=True)
            subprocess.run(["gunzip", "pdb_seqres.txt.gz"], check=True)
            subprocess.run(["makeblastdb", "-in", "pdb_seqres.txt", "-dbtype", "prot",
                            "-title", "PDBdb", "-parse_seqids", "-out", "PDBdb"], check=True)
            os.makedirs("db/blastp", exist_ok=True)

            files = glob.glob("PDBdb*")
            for filename in files:
                file_name = os.path.basename(filename)
                s = shutil.move(filename, f"db/blastp/{file_name}")

        # BLASTP
        blastp_cline = NcbiblastpCommandline(query = input_file, db = "db/blastp/PDBdb", outfmt = 5, out = "results.xml")
        stdout, stderr = blastp_cline()

        # Parse the results file
        E_VALUE_THRESH = 1e-20
        hits_dict = {}
        for record in NCBIXML.parse(open("results.xml", "r")):
            if record.alignments:
                for align in record.alignments:
                    for hsp in align.hsps:
                        if hsp.expect < E_VALUE_THRESH:
                            hits_dict.setdefault(align.title[4:8], align.title[9])
        protein_codes = list(hits_dict.keys())
        protein_codes = [x.lower() for x in protein_codes]
        protein_chains = list(hits_dict.values())

        if len(hits_dict) < 12: # arbitrary value to ensure we have an enough number of proteins
            print("Less than 7 proteins were found with that threshold, remote homologs will be searched... Introduce a higher threshold if only close homologs are desired")
            # Psiblast db
            if os.path.isdir('./db/psiblast') and os.path.isdir('./db/blastp'):
                if os.listdir("./db/psiblast") and os.listdir("./db/blastp"):
                    print("Using the existing PDB and psiblast databases located at './db/blastp' and './db/psiblast'")
                else:
                    raise OSError("Not database found! Please delete the ./db/blastp folder to allow the creation of a new one")
            else:
                print("Creating a Uniprot database to perform the psiblast")
                subprocess.run(["wget", "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"], check=True)
                subprocess.run(["gunzip", "uniprot_sprot.fasta.gz"], check=True)
                subprocess.run(["makeblastdb", "-in", "uniprot_sprot.fasta", "-dbtype", "prot",
                                "-title", "UNIPROTdb", "-parse_seqids", "-out", "UNIPROTdb"], check=True)
                os.makedirs("db/psiblast", exist_ok=True)

                files = glob.glob("UNIPROTdb*")
                for filename in files:
                    file_name = os.path.basename(filename)
                    s = shutil.move(filename, f"db/psiblast/{file_name}")

            # PSIBLAST
            psiblast_uniprot_cline = NcbipsiblastCommandline(
                db = 'db/psiblast/UNIPROTdb', query = input_file, evalue =  1 ,
                out = "psiblast_uniprot_results.xml", outfmt = 5,
                out_pssm = "psiblast_uniprot3.pssm", num_iterations = 3
                )
            stdout_psi_uni, stderr_psi_uni = psiblast_uniprot_cline()
            psiblast_pdb_cline = NcbipsiblastCommandline(
                db = 'db/blastp/PDBdb', out = "psiblast_pdb_results.xml", outfmt = 5,
                in_pssm = "psiblast_uniprot3.pssm"
                )
            stdout_psi_pdb, stderr_psi_pdb = psiblast_pdb_cline()
            PSI_E_VALUE_THRESH = 1e-4
            for psirecord in NCBIXML.parse(open("psiblast_pdb_results.xml", "r")):
                if psirecord.alignments:
                    for psialign in psirecord.alignments:
                        for hsp in psialign.hsps:
                            if hsp.expect < PSI_E_VALUE_THRESH:
                                hits_dict.setdefault(psialign.title[4:8], psialign.title[9])
                if len(hits_dict) == 12:
                    break
            protein_codes = list(hits_dict.keys())
            protein_codes = [x.lower() for x in protein_codes]
            protein_chains = list(hits_dict.values())

            return protein_codes, protein_chains, hits_dict

        else:
            return protein_codes, protein_chains, hits_dict

    else:
        print("Sorry, introduce a valid input")


def get_pdb_from_alphafold(target_id):
    """Retrieve structure from alphafold server
    Return filename
    """
    import requests

    url = f"https://alphafold.ebi.ac.uk/files/AF-{target_id}-F1-model_v2.pdb"
    # Define the local filename to save data
    local_file = f'structures/{target_id}.pdb'
    # Make http request for remote file data
    target_data = requests.get(url)
    # Save file data to local copy
    with open(local_file, 'wb')as file:
        file.write(target_data.content)

    return local_file

def get_pdb_structure(pdb_file, pdb_code):
    """ Parse pdb structure
    Return: structure object
    """
    from Bio.PDB import PDBParser

    # Read pdb file
    with open(pdb_file, "r") as fh:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_code, fh)

    return structure

def extract_fasta_sequence(fasta_file):
    """Parse FASTA file

    Return: protein id and protein sequence
    """
    from Bio import SeqIO

    with open(fasta_file, 'r') as fasta:
        for i, record in enumerate(SeqIO.parse(fasta, 'fasta')):
            if i == 1:
                raise OSError("The FASTA input file must contain one unique record")

            prot_id = record.id
            prot_seq = record.seq

    return prot_id, prot_seq

def pdb_quality_filter(protein_codes, hits_dict, resolution = 2):
    # Needs the protein_codes and hits_dict from the get_pdb_homologs function
    """
    Downloads a list of PDB files and stores them in the structures directory. Then, takes the files
    from the "structures/" directory, parses them to obtain the resolution and returns the three
    hits with the highest e-value and a resolution higher than 2 Angstroms, by default.

    Input: List of protein codes from the BLAST and dictionary of the chains for each hit. Optional: Resolution

    Return: List of three protein PDB codes and their respective chains. Directory called "structures" with
    the PDB files in it in the format: "structures/pdb{code}.ent"
    """
    # SHOULD RETURN THE OUTPUT FILES?

    from Bio.PDB.PDBList import PDBList
    import os

    final_protein_codes = []
    final_chains = []
    #pdb_outfiles = []
    for code in protein_codes:
        ## Get PDB files
        r = PDBList()  # Object instance
        r.retrieve_pdb_file(pdb_code = code, file_format = "pdb", pdir = "structures/", overwrite=True) # creates the directory if do not exists
        # Read pdb file
        structure = get_pdb_structure(f"structures/pdb{code}.ent", code)
        # Evaluate resolution
        if structure.header['resolution']:
            if structure.header['resolution'] <= resolution:
                final_protein_codes.append(code)
                final_chains.append(hits_dict[code])
                #pdb_outfiles.append(f"structures/pdb{code}.pdb") and include variable in return
            else:
                os.remove(f"structures/pdb{code}.ent")
        else:
            os.remove(f"structures/pdb{code}.ent")

        if len(final_protein_codes) == 3:
            break
    return final_protein_codes, final_chains

def get_aln_sequences(input_file, pdb_outfiles, pdb_codes, chains):
    """Construct a FASTA file containing the sequences of the target protein and its pdb homologs

    Return: Fasta filename
    """

    from Bio import SeqIO
    ## To track events
    import logging
    ## To hide warnings
    import warnings

    # 1. Obtain PDB files for candidates and target
    logging.info("Obtaining PDB files for candidate and target proteins")


    # Save the fasta_sequences
    fasta_outfile = "structures/pdb_sequences.fa"

    with warnings.catch_warnings(record=True) as caught_warnings:
        with open(fasta_outfile, "w") as out:
            # attach target sequence
            prot_id, prot_seq = extract_fasta_sequence(input_file)
            print(">"+prot_id, file = out)
            print(prot_seq, file = out)

            # attach sequence from pdb homologs
            for i in range(len(pdb_outfiles)):
                with open(pdb_outfiles[i], 'r') as pdb_fh:
                    for record in SeqIO.parse(pdb_fh, 'pdb-atom'):
                        if f":{chains[i]}" in record.id:
                            fasta_id = f">{record.id}"
                            fasta_seq = record.seq

                            print(fasta_id, file = out)
                            print(fasta_seq, file = out)

        return fasta_outfile

def msa(pdb_seqs_filename):
    """Perform multiple sequence alignment with Tcoffe

    Return: multiple Aligment object
    """

    from Bio.Align.Applications import TCoffeeCommandline
    from Bio import AlignIO

    # Perform  multiple sequence alignment (MSA) with T-Coffe
    msa_outfile = "structures/alignment.aln"
    tcoffee_cline = TCoffeeCommandline(infile=pdb_seqs_filename,
                                    output="clustalw",
                                    outfile=msa_outfile)
    tcoffee_cline()

    # Read the MSA
    msa_obj = AlignIO.read(msa_outfile, "clustal")

    return msa_obj

def get_modified_bfactors(msa_seqs, all_norm_bfactors, target_norm_bfactors, all_matrices):
    from statistics import mean

    flex_scores = {
        'A': 0.717,'C': 0.668,'D': 0.921,'E': 0.963,'F': 0.599, 'G': 0.843,
        'H': 0.754,'I': 0.632,'K': 0.912,'L': 0.681,'M': 0.685,'N': 0.851,
        'P': 0.850,'Q': 0.849,'R': 0.814,'S': 0.840,'T': 0.758,'V': 0.619,
        'W': 0.626,'Y': 0.615
        }

    # Target msa
    msa_seq_target = msa_seqs.pop(0)
    matrix_target = all_matrices.pop(0)

    # Obtain length of aln
    msa_len = len(msa_seq_target)

    # Number of homologues
    n_homologs = len(msa_seqs)

    # PDB residue index
    len_pdb = -1
    # PDB sequence
    seq = ""

    for i in range(msa_len):
        target_aa = msa_seq_target[i]
        # Missing residues (X) or gaps (-) are ignored
        if target_aa != "X" and target_aa != "-":
            len_pdb += 1
            # Check the same position in the homolog sequence
            h_aa =[homolog_seq[i] for homolog_seq in msa_seqs]

            if target_aa in h_aa:
                # If there are coincidences between target and homologues
                identities = h_aa.count(target_aa)

                if identities == 1:
                    # Find in which homologue the coincidence is found
                    homolog_index = h_aa.index(target_aa)

                    # Get position of this coincident aa in the homologue sequence and obtain the bfactor
                    residue_index = all_matrices[homolog_index][i]
                    bfactor = all_norm_bfactors[homolog_index][residue_index]
                else:
                    # If there are more than one coincidence
                    # Find in which homologue the coincidences are found
                    homolog_indexes = [i for i in range(len(h_aa)) if h_aa[i] == target_aa]

                    # For each homologue, find the position where the coincident is located in each homologue
                    residues_indexes = [all_matrices[homolog_index][i] for homolog_index in homolog_indexes]

                    # Get the bfactors for each coincidence
                    bfactors = [all_norm_bfactors[val[0]][val[1]] for val in list(zip(homolog_indexes, residues_indexes))]

                    #Compute the mean of b_factors
                    bfactor =  mean(bfactors)

            elif len(set(h_aa)) > 1 :
                h_aa_set = set(h_aa)

                # Remove possible X
                if 'X' in h_aa_set:
                    h_aa_set.remove('X')

                # If there there is some aa aligned with target
                no_aln_set = set('-')

                # Find the different residue
                diff_set = set(h_aa).difference(no_aln_set)

                if len(diff_set) == 1:
                    # Find in which homologue the residue is found
                    res = diff_set.pop()
                    homolog_index = h_aa.index(res)

                    # Get position of this coincident aa in the homologue sequence and obtain the bfactor
                    residue_index = all_matrices[homolog_index][i]
                    bfactor = all_norm_bfactors[homolog_index][residue_index]
                else:
                    # There are more than one different residue aligned with our target aa
                    # Use predefined bfactor values
                    bfactor = flex_scores[target_aa]

            else:
                # There is no homolog aa aligned with our target aa
                # Use predefined bfactor values
                bfactor = flex_scores[target_aa]

            target_norm_bfactors[len_pdb] = bfactor

    return target_norm_bfactors

def get_pdb_seq(msa_obj):
    """
    Return: pdb sequence len and pdb sequence
    """
    # Obtain length of aln
    msa_len = msa_obj.get_alignment_length()

    # Outputs
    all_seqs    = list()
    all_len_pdb = list()
    all_matrices = list()

    for obj in msa_obj:
        # Obtain  sequence
        msa_seq = obj.seq
        # PDB residue index
        len_pdb = -1
        seq = ''
        matrix = list()

        for i in range(msa_len):
            # Missing residues (X) or gaps (-) are ignored
            if msa_seq[i] != "X" and msa_seq[i] != "-":
                seq += msa_seq[i]
                len_pdb += 1

                matrix.append(len_pdb)
            else:
                matrix.append(None)

        all_seqs.append(seq)
        all_len_pdb.append(len_pdb)
        all_matrices.append(matrix)

    return all_seqs, all_len_pdb, all_matrices


def get_bfactors(pdb_file, pdb_code, pdb_chain, len_pdb):
    """Obtain the b-factors of c-alpha atoms for all the residues involved in the alignment

    Return: list of b-factors
    """

    # Read pdb file
    structure = get_pdb_structure(pdb_file, pdb_code)

    for chain in structure.get_chains():
        if chain.id == pdb_chain:
            chain_struct = chain
    # Get the residues involved in the alignment
    residues = [res for res in chain_struct.get_residues()]

    bfactors = list()
    hetatoms = 0
    for i in range(len(residues)):
        if i <= len_pdb:
            try:
                bfactors.append(residues[i].child_dict['CA'].get_bfactor())
            except KeyError:
                print(f"residue({i}) {residues[i]} has not CA. Probably an heteroatom found. Skipping..")
                hetatoms += 1
                continue
    if hetatoms:
        for i in range(len_pdb-1,len_pdb-1+hetatoms):
            bfactors.append(residues[i].child_dict['CA'].get_bfactor())

    return bfactors




def normalize_bfactor(bfactor):
    """Obtain the normalized b-factors of c-alpha atoms using the formula:
        Bnorm = (B-mean(B))/stedev(B)

    Return: list of normalized b-factors
    """
    import statistics
    n_bfactor = list(map(lambda bf:
                (bf-statistics.mean(bfactor))/statistics.stdev(bfactor), bfactor))
    return n_bfactor


## Flexibility
    # [?] Si tenemos b-values de nuestra proteina problema, yo aplicaría la fórmula
    # directamente
    # Si no los tenemos, he pensado en dos posibles approachs:
    # 1 APPROACH
    # Sobre lo que tenemos hecho hasta ahora, calcularia la media de los diferentes
    # valores del b-factor de un aa, dados por las proteinas homologas. Luego,
    # aplicaría la formula del paper https://www.sciencedirect.com/science/article/pii/S000527361300206X?via%3Dihub#bb0145
    # del apartado 2.3.
    # 2 APPROACH
    # Pillar un homólogo de referencia directamente y usar los valores de flexibilada
    # obtenidos en el paper https://onlinelibrary.wiley.com/doi/full/10.1110/ps.0236203
    # (los resultados están en las Tablas 3 y 4, de hecho en el primer paper
    # aplican la formula tomando como valores lambda los recogidos en este ultimo)

# FIRST APPROACH
def flexibility(final_bfactors, window_size):
    # [?] Valorar si incluir el scale_factor o no
    """ This function will return a flexibility score for the total protein and a
        list with the flexibility of each aa, taking into account how flexible or
        rigid are the neighbourhoods of each amino acid.
        We will define the flexibility index as a weighted average of amino acid
        flexibility over wole residue chain of the protein. The neighbourhoods effect
        is considered using a sliding hat-shaped window. The following formula is
        implemented:
         F_ws = scale_factor * sum(j=s,j=(L-(s-1))){lambda_j+sum(i=1,i=s-1){i/(s)*(lambda_{j-1}+lambda_{j+1})}}
        where:
        \lambda_j: mean of b-values for each set of 3(?) aminoacids
        s = (ws+1)/2 ("start index")
        L = length of the peptide (in this case aa of similarity regions)
        ws = window_size
        scale_factor = 1/(s(L-(ws-1)))

        OUTPUT: protein flexibility score, list with the flexibility score of each aa
    """
    import statistics

    # Length of the peptidic chain
    L = len(final_bfactors)
    s = int((window_size + 1)/2)
    scale_factor = 1/s*(L-(window_size-1))
    # Compute the FORMULA
    all_aa_flex = []
    for j in range(s-1,L-(s-1)):
        k = 0
        for i in range(1,s-1):
            k = i/s*(final_bfactors[j-1]+final_bfactors[j+1])
        #Flexibility score F of each aa
        all_aa_flex.append(final_bfactors[j] + k)


    # Flexibility score F of the whole protein
    # With scale factor
    F = scale_factor*sum(all_aa_flex)
    # Without scale factor
    # F = sum(all_aa_flex)
    return F, all_aa_flex

def scale_function(all_aa_flex):
    """Obtain the flexibility aminoacids scores in a range of values between 0 and 1,
    following the formula:

        Fscore = (Fscore-min(Fscore))/(max(Fscore)-min(Fscore))

    Return: list of scaled F scores
    """
    all_aa_flex_norm = list(map(lambda i:(i - min(all_aa_flex))/(max(all_aa_flex) - min(all_aa_flex)), all_aa_flex))
    return all_aa_flex_norm

def get_sstructure(pdb_file, pdb_code, prot_chain=''):
    """Calculate the secondary structure and accesibility by using DSSP program
        - Needs pdb extension


        The DSSP codes for secondary structure used here are:
            =====     ====
            Code      Structure
            =====     ====
             H        Alpha helix (4-12)
             B        Isolated beta-bridge residue
             E        Strand
             G        3-10 helix
             I        Pi helix
             T        Turn
             S        Bend
             \-       None
            =====     ====

        However, we will convert DSSP's 8-state assignments into 3-state:
        [C - coil, E - extended (beta-strand), H - helix].

    Return: string with secondary structure
    """
    from Bio.PDB.DSSP import DSSP

    # Parse PDB file and execute DSSP
    structure = get_pdb_structure(pdb_file, pdb_code)
    model = structure[0]
    dssp = DSSP(model, pdb_file,  dssp='mkdssp')

    # Retrieve the secondary structure from dssp results
    sstructure = ''
    for i in range(len(dssp)):
        a_key = list(dssp.keys())[i]

        if prot_chain:
            if a_key[0] == prot_chain:
                sstructure += dssp[a_key][2]
        else:
            sstructure += dssp[a_key][2]

    # Convert 8-state to 3-state
    sstructure = sstructure.replace('-', 'C')
    sstructure = sstructure.replace('I', 'C')
    sstructure = sstructure.replace('T', 'C')
    sstructure = sstructure.replace('S', 'C')
    sstructure = sstructure.replace('G', 'H')
    sstructure = sstructure.replace('B', 'E')

    return sstructure

def get_hydrophobicity(fasta_file):
    """
    Calculates the hydrophobicity score for each aminoacid and for the global protein (GRAVY)

    Input: FASTA file with the sequence

    Return: List of scores per aminoacid and GRAVY score
    """
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    from Bio.SeqUtils.ProtParam import ProtParamData

    # Parse Fasta File
    prot_id, prot_seq = extract_fasta_sequence(fasta_file)
    prot_seq = str(prot_seq)

    # Obtain hydrophobicity scores
    X = ProteinAnalysis(prot_seq)
    hydroph_scores_aa = X.protein_scale(window=7, param_dict=ProtParamData.kd)
    gravy = sum(hydroph_scores_aa)/len(prot_seq)

    return hydroph_scores_aa, gravy

def from_sstructure_to_score(string_sstructure):
    """
    Transforms the string of secondary structures to scores for the flexibility plot.
    Input: String with 3-state secondary structures
    Return: List of scores for the secondary structure
    """
    list_sstructure = list(string_sstructure)
    for index, item in enumerate(list_sstructure):
        if item == "H":
            list_sstructure[index] = 0
        elif item == "E":
            list_sstructure[index] = 0.5
        elif item == "C":
            list_sstructure[index] = 1
    return list_sstructure

def data_frame_results(norm_flex_scores,hydroph_scores, target_seq, list_sstructure, sstructures, ws):
    """ Functions to save the amino acids sequence, flexibility score, hydrophobicity
    score and secondary structure of the protein in a DataFrame. The function takes
    into account the limitation of losing amino acids with the window_size.
    It will return two DataFrames, one to write the output file (df_results)
    with the "origianl" secondary structure, and the other one to use in the plot
    function with the secondary structure transformed to numberes (df_plot).

    Output: DataFrames df_results and df_plot with results.
    """
    import pandas as pd
    # 1. Save the needed values info
    i = int(( ws -1 )/2)
    L = len(target_seq)
    amino_acids = list(target_seq)[i:L-i]
    sstructure = list_sstructure[i:L-i]
    hydroph_scores_reverse = [num * -1 for num in hydroph_scores]

    # 2. Create the DataFrames
    ## DataFrame adapted to plot the results:
    df_plot = pd.DataFrame(list(zip(sstructure, hydroph_scores_reverse, norm_flex_scores, amino_acids)),
    columns = ["sstructure","hydrophobicity", "flex_scores", "amino_acids"])
    ## DataFrame with "real" results
    df_results = pd.DataFrame(list(zip(amino_acids, norm_flex_scores, hydroph_scores, sstructures)),
    columns = ["amino_acids","flex_scores","hydrophobicity","sstructure"])

    return df_results, df_plot

def plot_heatmap(ax, cmap, col, df_short, i, L = 0):
    """ Function that plots the flexibility scores, the hydrophobicity scores and
    the secondary structure by amino acid. Darker colors are asociated with less
    flexible and less hydrophobic zones. In the same way, darker color is assigned
    to helix structure and lighter one to coil structure.
    INPUT:
            ax: argument that will represent the plot
            cmap: scale_colors
            col: data frame column names with the names features to represent
            df_short: df contianing part of the whole df (50 values in this case)
            i: index that represents the splitted df that is being plotted
            L: len of the whole data frame

    OUTPUT: Plot with flexibility scores, hydrophobicity values and secondary structure.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    # 1. Save information from the df_short
    aa = df_short["amino_acids"]
    l = len(df_short)

    # 2. Create the plot:
    ## Iterate over features
    for j in range(0,3):

        x = list(range(0,l))
        y = [j]*len(x)
        color = cmap(df_short[col[j]])
        p = ax.scatter(x, y, color = color, cmap=cmap, s = 120) # s = shape

    ## Remove all spines
    ax.set_frame_on(False)
    ## Set grid lines with some transparency
    ax.grid(alpha=0.4)
    ## Make sure grid lines are behind other objects
    ax.set_axisbelow(True)
    ## Set position for x ticks
    ax.set_xticks(np.arange(len(aa)))
    ## Set labels for the x ticks (the names of the types of plastic)
    ax.set_xticklabels(aa)
    ## Set position for y ticks
    ax.set_yticks(np.arange(len(col[:3])))
    ## Set labels for the y ticks (the names of the types of plastic)
    ax.set_yticklabels(col[:3])
    ## Remove tick marks by setting their size to 0. Set text color to "0.3" (a type of grey)
    ax.tick_params(size=0, colors="0.3")
    ## Set label for horizontal axis.
    ax.set_xlabel("Analysis", loc="center")
    ## Add leyend

    return p

def plot_linear(df_plot):
    """
    """
    import seaborn as sns
    sns.set_theme(style="darkgrid")
    L = len(df_plot)
    df_plot.insert(4, "aa_residue_num", range(1,L+1))
    p = sns.lineplot(x = "aa_residue_num", y = "flex_scores",
             data=df_plot)
    p.set_title("Flex_scores distribution")
