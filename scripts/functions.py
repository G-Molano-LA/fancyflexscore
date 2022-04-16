def get_pdb_homologs(input_file, E_VALUE_THRESH = 1e-20):
    """Obtain PDB codes and chains from protein homologs of the query protein

    This function conducts a BLASTP against the PDB database locally and retrieves the first 3 unique
    hits with their respective chains. Proteins with homology in several chains are not considered unique
    and only the first chain is kept.
    In case less than 3 hits, with a lower threshold than 1e-20, are found two psi-blasts are conducted. First,
    against the Uniprot database to obtain the 3rd iteration PSSM and then against the PDB database using the
    obtained PSSM. This is only conducted to fill the three homologs in case not enough homologs are obtained
    with the BLASTP.

    Input: FASTA file with the query protein sequence

    Return: List of protein codes and list of protein chains of the three homologs.
    [ ] Maybe, we can return the whole list of hits, in order to be able to access to it later
    """

    # Python modules
    import os
    from Bio import SeqIO
    from Bio.Blast.Applications import NcbiblastpCommandline
    from Bio.Blast import NCBIXML
    from Bio.Blast.Applications import NcbipsiblastCommandline

    # Checking if the input is a file or a sequence:
    if os.path.isfile(input_file):
        # BLASTP
        blastp_cline = NcbiblastpCommandline(query = input_file, db = "../db/PDBdb", outfmt = 5, out = "results.xml")
        stdout, stderr = blastp_cline()

        # Parse the results file
        # [?] Poner e-value threshold como argumento modificable del usuario
        # [?] Porque esta en mayúscula?
        E_VALUE_THRESH = 1e-20
        hits_dict = {}
        for record in NCBIXML.parse(open("results.xml", "r")):
            if record.alignments:
                for align in record.alignments:
                    for hsp in align.hsps:
                        if hsp.expect < E_VALUE_THRESH:
                            hits_dict.setdefault(align.title[4:8], align.title[9])

        # Obtain just the first three unique hits
        if len(hits_dict) >= 3:
            protein_codes = list(hits_dict.keys())[:3]
            protein_codes = [x.lower() for x in protein_codes]
            protein_chains = list(hits_dict.values())[:3]
        else:
            # PSIBLAST
            # [?] Poner e-value threshold como argumento modificable del usuario
            # [?] Porque esta en mayúscula?
            psiblast_uniprot_cline = NcbipsiblastCommandline(
                db = '../db/UNIPROTdb', query = input_file, evalue =  1 ,
                out = "psiblast_uniprot_results.xml", outfmt = 5,
                out_pssm = "psiblast_uniprot3.pssm", num_iterations = 3
                )
            stdout_psi_uni, stderr_psi_uni = psiblast_uniprot_cline()
            psiblast_pdb_cline = NcbipsiblastCommandline(
                db = '../db/PDBdb', out = "psiblast_pdb_results.xml", outfmt = 5,
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
                if len(hits_dict) == 3:
                    break
            protein_codes = list(hits_dict.keys())[:3]
            protein_codes = [x.lower() for x in protein_codes]
            protein_chains = list(hits_dict.values())[:3]

        return protein_codes, protein_chains

    else:
        print("Sorry, introduce a valid input")


def get_pdb_sequences(pdb_codes, chains, pdb_outfiles):
    """Obtain PDB files and fasta sequences for protein homologues

    The pdb files and fasta sequences are stores in 'structures' directory.
    Generates a FASTA file containningg all the sequences of the pdb structures.

    PDB generated files: f"structures/pdb{code}.ent"

    Return: Fasta filename
    """

    from Bio.PDB.PDBList import PDBList   # to get PDB file from PDB server
    from Bio import SeqIO
    ## To track events
    import logging

    # 1. Obtain PDB files for candidates and target
    logging.info("Obtaining PDB files for candidate and target proteins")

    ## Get PDBs
    r = PDBList()  # Object instance
    r.download_pdb_files(pdb_codes = pdb_codes , file_format = "pdb",
                         pdir = "structures/", overwrite=True)   # creates the directory if do not exists

    # Save the pdb sequences into a file
    outfile = "structures/pdb_sequences.fa"
    with open(outfile, "w") as out:
        for i in range(len(pdb_outfiles)):
            with open(pdb_outfiles[i], 'r') as pdb_fh:
                for record in SeqIO.parse(pdb_fh, 'pdb-atom'):
                    if f":{chains[i]}" in record.id:
                        fasta_id = f">{record.id}"
                        fasta_seq = record.seq

                        print(fasta_id, file = out)
                        print(fasta_seq, file = out)

        return outfile

def msa(pdb_seqs_filename):
    """Perform multiple sequence alignment with Tcoffe
    [ ] Include the input file

    Return: multiple Aligment object
    """

    # install t-coffe: https://www.tcoffee.org/Projects/tcoffee/documentation/index.html#document-tcoffee_installation

    from Bio.Align.Applications import TCoffeeCommandline
    from Bio import AlignIO

    # Perform  multiple sequence alignment (MSA) with T-Coffe
    msa_outfile = "structures/alignment.aln"
    tcoffee_cline = TCoffeeCommandline(infile=pdb_seqs_filename,
                                    output="clustalw",
                                    outfile=msa_outfile)
    tcoffee_cline()

    # Read the MSA
    msa_obj             = AlignIO.read(msa_outfile, "clustal")

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
    # all then position of similarities for each protein
    all_sim = list()
    all_aln_res = list()

    for obj in msa_obj:
        # Obtain  sequence
        msa_seq = obj.seq
        # PDB residue index
        aln_res = -1
        # similarities for the sequence
        sim = list()

        for i in range(msa_len):
            # Missing residues (X) or gaps (-) are ignored
            if msa_seq[i] != "X" and msa_seq[i] != "-":
                aln_res += 1
            if annot[i] == "*":
                sim.append(aln_res)

        all_sim.append(sim)
        all_aln_res.append(aln_res)

    return all_sim, all_aln_res


def get_bfactors(pdb_file, pdb_code, pdb_chain, aln_res):
    """Obtain the b-factors of c-alpha atoms for all the residues involved in the alignment

    Return: list of b-factors
    """
    from Bio.PDB import PDBParser

    # Read pdb file
    with open(pdb_file, "r") as fh:
        parser = PDBParser()
        structure = parser.get_structure(pdb_code, fh)

    resolution = structure.header['resolution']

    for chain in structure.get_chains():
        if chain.id == pdb_chain:
            chain_struct = chain
    # Get the residues involved in the alignment
    residues = [res for res in chain_struct.get_residues()]

    # Get the bfactors of the pdb residues involved in the alignment
    bfactors = [residues[i].child_dict['CA'].get_bfactor() for i in range(len(residues))
                if i <= aln_res]
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
def flexibility(all_sim_bfactors, window_size):
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
    # Compute the mean of b_values
    mean_bvalues = [statistics.mean(i) for i in list(zip(all_sim_bfactors[0],all_sim_bfactors[1],all_sim_bfactors[2]))]
    # Length of the peptidic chain
    L = len(mean_bvalues)
    s = int((window_size + 1)/2)
    scale_factor = 1/s*(L-(window_size-1))
    # Compute the FORMULA
    all_aa_flex = []
    for j in range(s,L-(s-1)):
        k = 0
        for i in range(1,s-1):
            k = i/s*(mean_bvalues[j-1]+mean_bvalues[j+1])
        #Flexibility score F of each aa
        all_aa_flex.append(mean_bvalues[j-1] + k)


    # Flexibility score F of the whole protein
    # With scale factor
    F = scale_factor*sum(all_aa_flex)
    # Without scale factor
    # F = sum(all_aa_flex)
    return F, all_aa_flex

## SECOND APPROACH

def get_bfactors2(pdb_file, pdb_code, pdb_chain, all_sim, aln_res):
    """Obtain the b-factors of c-alpha atoms for all the residues involved in the alignment

    Return: list of b-factors
    """
    from Bio.PDB import PDBParser

    # Read pdb file
    with open(pdb_file, "r") as fh:
        parser = PDBParser()
        structure = parser.get_structure(pdb_code, fh)

    resolution = structure.header['resolution']

    for chain in structure.get_chains():
        if chain.id == pdb_chain:
            chain_struct = chain
    # Get the residues involved in the alignment
    residues = [res for res in chain_struct.get_residues()]

    # Get the residues of the pdb residues involved in the alignment
    residues_aln = [residues[i] for i in all_sim
                if i <= aln_res]

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


if __name__ == '__main__':
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
    ## Get the b-factors for all residues and normalize

    all_bfactors = [get_bfactors(prot[0], prot[1], prot[2], prot[3]) for prot in
                     list(zip(pdb_outfiles,protein_codes, protein_chains, all_aln_res))]
    all_norm_bfactors = [normalize_bfactor(bfactors) for bfactors in all_bfactors]

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
