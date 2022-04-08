
def get_pdb_homologs(input_file):
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
    """

    # Python modules
    import os
    from Bio import SeqIO
    from Bio.Blast.Applications import NcbiblastpCommandline
    from Bio.Blast import NCBIXML
    from Bio.Blast.Applications import NcbipsiblastCommandline
    # Checking if the input is a file or a sequence:
    if input_file:
        if os.path.isfile(input_file):
            blastp_cline = NcbiblastpCommandline(query = input_file, db = "../db/PDBdb", outfmt = 5, out = "results.xml")
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

            # Obtain just the first three unique hits
            if len(hits_dict) >= 3:
                protein_codes = list(hits_dict.keys())[:3]
                protein_codes = [x.lower() for x in protein_codes]
                protein_chains = list(hits_dict.values())[:3]
            else:
                psiblast_uniprot_cline = NcbipsiblastCommandline(db = '../db/UNIPROTdb', query = input_file, evalue =  1 , out = "psiblast_uniprot_results.xml", outfmt = 5, out_pssm = "psiblast_uniprot3.pssm", num_iterations = 3)
                stdout_psi_uni, stderr_psi_uni = psiblast_uniprot_cline()
                psiblast_pdb_cline = NcbipsiblastCommandline(db = '../db/PDBdb', out = "psiblast_pdb_results.xml", outfmt = 5, in_pssm = "psiblast_uniprot3.pssm")
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
    r = PDBList() # Object instance
    r.download_pdb_files(pdb_codes = pdb_codes , file_format = "pdb",
                        pdir = "structures/", overwrite=True) # creates the directory if do not exists

    # Save the pdb sequences into a file
    outfile = "structures/pdb_sequences.fa"
    with open(outfile, "w") as out:
        sequence_records = []
        for i in range(len(pdb_outfiles)):
            with open (pdb_outfiles[i], 'r') as pdb_fh:
                for record in SeqIO.parse(pdb_fh, 'pdb-atom'):
                    if f":{chains[i]}" in record.id:
                            fasta_id = f">{record.id}"
                            fasta_seq = record.seq

                            print(fasta_id, file = out)
                            print(fasta_seq, file = out)

        return outfile

def msa(pdb_seqs_filename):
    """Perform multiple sequence alignment with Tcoffe

    Return:Dictionary containing the position (key) and symbol (value) of aminoacid similarities
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
    msa             = AlignIO.read(msa_outfile, "clustal")
    # Obtain the sequence of reference
    msa_ref_seq     = msa[0].seq
    # Obtain the annotation of clustalw format
    msa_annotation  = read_clustal_annotations(msa_outfile)

    # Save the identifies in a list of tuples (aminoacid, position)
    ## Use 1-based position annotation (be consistent with pdb file)
    similarities = {"position": list(),
                    "symbol": list()}

    for i in range(len(msa_annotation)):
        if msa_annotation[i] == "*":
            similarities["position"].append(i+1)
            similarities["symbol"].append(msa_ref_seq[i])

    return similarities

def read_clustal_annotations(msa_outfile):
    """Obtain the annotations ('*', '.' , ':', ' ') of clustalw alignment format

    Return: annotations
    """
    with open(msa_outfile, 'r') as fh:
        annotation_line = ""

        for line in fh:
            if line.isspace():
                pass
            if '*' in line:
                line = line.strip()
                annotation_line += line
        return annotation_line

def get_bfactors(pdb_file, pdb_code, similarities):
    """Obtain the b-factors of c-alpha atoms for regions of similarities

    Return: list of b-factors
    """
    from Bio.PDB import PDBParser

    with open(pdb_file, "r") as fh:
        parser = PDBParser()
        structure = parser.get_structure(pdb_code, fh)

    bfactor = []
    for atom in structure.get_atoms():
        if atom.get_name() == 'CA':
            position = atom.get_serial_number()
            if position in similarities["position"]:
                bfactor.append(atom.get_bfactor())
    return bfactor

if __name__ == '__main__':
    protein_codes = ["1a3n", "2dhb"] # PDB codes in lowercasese
    protein_chains = ["A", "A"]
    # Output format of pdb_files
    pdb_outfiles = [f"structures/pdb{code}.ent" for code in protein_codes]


    pdb_seqs_filename = get_pdb_sequences(protein_codes, protein_chains, pdb_outfiles)
    similarities= msa(pdb_seqs_filename)

    bfactor = [get_bfactors(pdb_file, pdb_code, similarities) for pdb_file,pdb_code in zip(pdb_outfiles, protein_codes)]
