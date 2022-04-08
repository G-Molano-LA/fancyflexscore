
def msa(protein_codes, protein_chains):

    # install t-coffe: https://www.tcoffee.org/Projects/tcoffee/documentation/index.html#document-tcoffee_installation

    from Bio.Align.Applications import TCoffeeCommandline
    import get_PDB

    pdb_seqs_filename = get_PDB.get_pdb_sequences(protein_codes, protein_chains)

    tcoffee_cline = TCoffeeCommandline(infile=pdb_seqs_filename,
                                    output="clustalw",
                                    outfile="structures/alignment.aln")
    tcoffee_cline()


if __name__ == '__main__':
    protein_codes = ["1a3n", "2dhb"] # PDB codes in lowercasese
    protein_chains = ["A", "A"]
    msa(protein_codes, protein_chains)
