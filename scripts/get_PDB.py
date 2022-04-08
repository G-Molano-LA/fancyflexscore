

def get_pdb_sequences(pdb_codes, chains):
    """Obtain PDB files and fasta sequences for protein homologues

    The pdb files and fasta sequences are stores in 'structures' directory.
    Generates a FASTA file containningg all the sequences of the pdb structures.

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
    # Output format of pdb_files
    pdb_files = [f"structures/pdb{code}.ent" for code in pdb_codes]

    # Save the pdb sequences into a file
    outfile = "structures/pdb_sequences.fa"
    with open(outfile, "w") as out:
        sequence_records = []
        for i in range(len(pdb_files)):
            with open (pdb_files[i], 'r') as pdb_fh:
                for record in SeqIO.parse(pdb_fh, 'pdb-atom'):
                    if f":{chains[i]}" in record.id:
                            fasta_id = f">{record.id}"
                            fasta_seq = record.seq

                            print(fasta_id, file = out)
                            print(fasta_seq, file = out)

        return outfile

def msa(protein_codes, protein_chains):

    # install t-coffe: https://www.tcoffee.org/Projects/tcoffee/documentation/index.html#document-tcoffee_installation

    from Bio.Align.Applications import TCoffeeCommandline
    from Bio import AlignIO

    pdb_seqs_filename = get_pdb_sequences(protein_codes, protein_chains)

    msa_outfile = "structures/alignment.aln"
    tcoffee_cline = TCoffeeCommandline(infile=pdb_seqs_filename,
                                    output="clustalw",
                                    outfile=msa_outfile)
    tcoffee_cline()

    msa = AlignIO.read(msa_outfile, "clustal")
    
    return msa

if __name__ == '__main__':
    protein_codes = ["1a3n", "2dhb"] # PDB codes in lowercasese
    protein_chains = ["A", "A"]
    align = msa(protein_codes, protein_chains)
