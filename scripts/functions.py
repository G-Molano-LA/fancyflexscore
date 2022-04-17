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

def download_pdb_files(pdb_codes_list):
    """Downloads a list of PDB files and stores them in the structures directory

    Input: List of PDB codes of the proteins

    Return: Directory called "structures" with the PDB files in it in the format: "structures/pdb{code}.ent"
    """
    from Bio.PDB.PDBList import PDBList
    ## Output files
    pdb_outfiles = [f"structures/pdb{code}.pdb" for code in pdb_codes_list]

    ## Get PDB files
    r = PDBList()  # Object instance
    r.download_pdb_files(pdb_codes = pdb_codes_list , file_format = "pdb",
                        pdir = "structures/", overwrite=True)   # creates the directory if do not exists

    return pdb_outfiles

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
        for i, record in enumerate(SeqIO.parse(fasta, 'pdb-atom')):
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
        structure = get_pdb_structure(f"structures/pdb{code}.pdb", code)
        # Evaluate resolution
        if structure.header['resolution'] <= 2:
            final_protein_codes.append(code)
            final_chains.append(hits_dict[code])
            #pdb_outfiles.append(f"structures/pdb{code}.pdb") and include variable in return
        else:
            os.remove(f"structures/pdb{code}.pdb")
        
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
            print(prot_id, file = out)
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

def clustal_annotations(msa_obj):
    """Obtain all the residue positions that are identifical in the alignment.
     Clustalw alignment annotations:
        '*'  -- all residues or nucleotides in that column are identical
        ':'  -- conserved substitutions have been observed
        '.'  -- semi-conserved substitutions have been observed
        ' '  -- no match.

    Return:
        - List of list, containing the positions of similarity for each structure
        - List containing the residue length (position?) associated to the pbd sequence used in the alignment
    """
    # Get annotations of clustalw aln
    annot = msa_obj.column_annotations['clustal_consensus']

    # Obtain length of aln
    msa_len = msa_obj.get_alignment_length()
    # Obtain sequence
    all_seqs = list()
    # all then position of similarities for each protein
    all_sim = list()
    all_con = list()
    all_rest = lit()
    all_len_pdb = list()

    for i in len(range(msa_obj)):

        # Obtain  sequence
        msa_seq = msa_obj[i].seq
        # PDB residue index
        len_pdb = -1
        # PDB sequence
        seq = ""
        # similarities for the sequence
        sim = list() # identities (*)
        con = list() # conserved (:)
        rest = list() # rest(., ' ')

        for i in range(msa_len):
            # Missing residues (X) or gaps (-) are ignored
            if msa_seq[i] != "X" and msa_seq[i] != "-":
                seq += msa_seq[i]
                len_pdb += 1

                if annot[i] == "*":
                    sim.append(len_pdb)
                elif annot[i] == ":":
                    con.append(len_pdb)
                else:
                    rest.append(len_pdb)

        all_sim.append(sim)
        all_con.append(con)
        all_rest.append(rest)
        all_len_pdb.append(len_pdb)
        all_seqs.append(seq)

    return all_sim, all_con, all_rest, all_len_pdb, all_seqs


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

    # Get the bfactors of the pdb residues involved in the alignment
    bfactors = [residues[i].child_dict['CA'].get_bfactor() for i in range(len(residues))
                if i <= len_pdb]
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

def get_modified_bfactors(all_norm_bfactors, all_sim, all_con, all_rest, all_seqs):
    """ Get the b-factors modified depending if they have identify, conserved regions or something else
    in the aln.
        - Identity : compute the mean of b-factors
        - Conserved: obtain the b-factor of coincident residue
        - Rest: put iria's bvalues

    """
    from statistics import mean

    # 1. Get the values for the target seq
    target_norm_bfactors = all_norm_bfactors.pop(0)
    target_seq       = all_seqs.pop(0)
    target_sim       = all_sim.pop(0)
    target_con       = all_con.pop(0)
    target_rest      = all_rest.pop(0)


    # 2. Get the b-factor only for regions of similarity
    all_sim_bfactors = list()

    ## Obtain the b-values of regions of similarity for each homologue
    for val in list(zip(all_norm_bfactors, all_sim)):
        norm_bfactors = val[0]
        similarities  = val[1]

        sim_bfactors = [norm_bfactors[i] for i in similarities]
        all_sim_bfactors.append(sim_bfactors)

    ## Compute the mean of b_values
    mean_sim_bvalues = [mean(i) for i in list(zip(*all_sim_bfactors))]

    ## Assign the newbvalues to the target
    for val in list(zip(target_sim, mean_sim_bvalues)):
        i       = val[0]
        bfactor = val[1]

        target_norm_bfactors[i] = bfactor

    # 3. Get the bfactor of coincident residue for conserved regions
    for val in list(zip(all_norm_bfactors, all_con, all_seqs)):
        homolog_norm_bfactors  = val[0]
        homolog_conserved      = val[1]
        homolog_seq            = val[2]

        for i in range(len(target_con)):
            t = target_con[i]
            h = homolog_conserved[i]

            if target_seq[t] == homolog_seq[h]:
                target_norm_bfactors[t] = homolog_norm_bfactors[h]

    # 4. Get the bfactors of rest residues for non-conserved regions
    flex_scores = {
        'A': 0.717,'C': 0.668,'D': 0.921,'E': 0.963,'F': 0.599, 'G': 0.843,
        'H': 0.754,'I': 0.632,'K': 0.912,'L': 0.681,'M': 0.685,'N': 0.851,
        'P': 0.850,'Q': 0.849,'R': 0.814,'S': 0.840,'T': 0.758,'V': 0.619,
        'W': 0.626,'Y': 0.615
        }

    for i in target_con:
        aa = target_seq[i]
        target_norm_bfactors[i] = flex_scores[aa]

    return target_norm_bfactors

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
    for j in range(s,L-(s-1)):
        k = 0
        for i in range(1,s-1):
            k = i/s*(final_bfactors[j-1]+final_bfactors[j+1])
        #Flexibility score F of each aa
        all_aa_flex.append(final_bfactors[j-1] + k)


    # Flexibility score F of the whole protein
    # With scale factor
    F = scale_factor*sum(all_aa_flex)
    # Without scale factor
    # F = sum(all_aa_flex)
    return F, all_aa_flex

## SECOND APPROACH

def get_bfactors2(pdb_file, pdb_code, pdb_chain, all_sim, len_pdb):
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

    # Get the residues of the pdb residues involved in the alignment
    residues_aln = [residues[i] for i in all_sim
                if i <= len_pdb]

    return residues_aln

def flexibility_val(res_names, window_size):
    #[?] Decide which flex score from the paper (Table 3 or 4) to use
    """ This function will return a flexibility score for the total protein and a
        list with the flexibility of each aa, taking into account how flexible or
        rigid are the neighbourhoods of each amino acid.
        We will define the flexibility index as a weighted average of amino acid
        flexibility over wole residue chain of the protein. The neighbourhoods effect
        is considered using a sliding hat-shaped window. The following formula is
        implemented:
         F_ws = scale_factor * sum(j=s,j=(L-(s-1))){lambda_j+sum(i=1,i=s-1){i/(s)*(lambda_{j-1}+lambda_{j+1})}}
        where:
        \lambda_j: flexibility scores obtained by David K. Smith,Predrag Radivojac,Zoran Obradovic,A. Keith Dunker,Guang Zhu
        s = (ws+1)/2 ("start index")
        L = length of the peptide (in this case aa of similarity regions)
        ws = window_size
        scale_factor = 1/(s(L-(ws-1)))

        OUTPUT: protein flexibility score, list with the flexibility score of each aa
    """

    # FLEXIBILITY scores
    flex_scores = {'ALA': 0.717,'CYS': 0.668,'ASP': 0.921,'GLU': 0.963,'PHE': 0.599,'GLY': 0.843,'HIS': 0.754,'ILE': 0.632,'LYS': 0.912,'LEU': 0.681,'MET': 0.685,'ASN': 0.851,'PRO': 0.850,'GLN': 0.849,'ARG': 0.814,'SER': 0.840,'THR': 0.758,'VAL': 0.619,'TRP': 0.626,'TYR': 0.615}
    flex_scores2 = {'ALA':-0.605,'CYS':-0.692,'ASP':-0.279,'GLU':-0.160,'PHE':-0.719,'GLY':-0.537,'HIS':-0.662,'ILE':-0.682,'LYS':-0.043,'LEU':-0.631,'MET':-0.626,'ASN':-0.381,'PRO':-0.271,'GLN':-0.368,'ARG':-0.448,'SER':-0.424,'THR':-0.525,'VAL':-0.669,'TRP':-0.727,'TYR':-0.721}
    # Length of the peptidic chain
    L = len(res_names)
    # Get the flex scores
    flex_scores_peptide = [flex_scores2[res_names[i]] for i in range(L)]
    s = int((window_size + 1)/2)
    scale_factor = 1/s*(L-(window_size-1))
    # Compute the FORMULA
    all_aa_flex_val = []
    for j in range(s,L-(s-1)):
        k = 0
        for i in range(1,s-1):
            k = i/s*(flex_scores_peptide[j-1]+flex_scores_peptide[j+1])
        #Flexibility score F of each aa
        all_aa_flex_val.append(flex_scores_peptide[j-1] + k)


    # Flexibility score F of the whole protein
    # With scale factor
    F_val = scale_factor*sum(all_aa_flex_val)
    # Without scale factor
    # F_val = sum(all_aa_flex)
    return F_val, all_aa_flex_val

def scale_function(all_aa_flex):
    """Obtain the flexibility aminoacids scores in a range of values between 0 and 1,
    following the formula:

        Fscore = (Fscore-min(Fscore))/(max(Fscore)-min(Fscore))

    Return: list of scaled F scores
    """
    all_aa_flex_norm = list(map(lambda i:(i - min(all_aa_flex))/(max(all_aa_flex) - min(all_aa_flex)), all_aa_flex))
    return all_aa_flex_norm

def get_sstructure(pdb_file, pdb_code, prot_chain):
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

    Return: sequence and sstructure
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

        if a_key[0] == prot_chain:
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

    # Obtain hydrophobicity scores
    X = ProteinAnalysis(prot_seq)
    hydroph_scores_aa = X.protein_scale(window=9, param_dict=ProtParamData.kd)
    gravy = sum(hydroph_scores_aa)/len(prot_seq)

    return hydroph_scores_aa, gravy

def from_sstructure_to_score(string_sstructure):
    """
    Transforms the string of secondary structures to scores for the flexibility plot.

    Input: String with 3-state secondary structures

    Return: List of scores for the secondary structure
    """
    string_sstructure = string_sstructure.replace('H','0')
    string_sstructure = string_sstructure.replace('E','1')
    string_sstructure = string_sstructure.replace('C','2')
    list_sstructure = list(string_sstructure)
    return list_sstructure
