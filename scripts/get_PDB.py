def get_pdb_sequences(pdb_codes, chains):
    """Obtain PDB files and fasta sequences for protein homologues

    The pdb files and fasta sequences are stores in 'structures' directory.

    Return: fasta sequences?
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


    for i in range(len(pdb_files)):
        with open (pdb_files[i], 'r') as pdb_fh:
            for record in SeqIO.parse(pdb_fh, 'pdb-atom'):
                if f":{protein_chains[i]}" in record.id:
                    fasta_id = f">{record.id}")
                    fasta_seq = record.seq)

if __name__ == '__main__':
    protein_codes = ["1a3n", "2dhb"] # PDB codes in lowercasese
    protein_chains = ["A", "A"]
    my_prots = get_pdb_sequences(protein_codes, protein_chains)
