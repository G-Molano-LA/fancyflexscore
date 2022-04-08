

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
