
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
    from Bio import SeqIO
    from Bio.Blast.Applications import NcbiblastpCommandline
    from Bio.Blast import NCBIXML
    from Bio.Blast.Applications import NcbipsiblastCommandline

    # Checking if the input is a file or a sequence:
    if os.path.isfile(input_file):
        # BLASTP
        blastp_cline = NcbiblastpCommandline(query = input_file, db = "db/PDBdb", outfmt = 5, out = "results.xml")
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
            
        if len(hits_dict) < 7: # arbitrary value to ensure we have an enough number of proteins
            print("Less than 7 proteins were found with that threshold, remote homologs will be searched... Introduce a higher threshold if only close homologs are desired")
            # PSIBLAST
            psiblast_uniprot_cline = NcbipsiblastCommandline(
                db = 'db/UNIPROTdb', query = input_file, evalue =  1 ,
                out = "psiblast_uniprot_results.xml", outfmt = 5,
                out_pssm = "psiblast_uniprot3.pssm", num_iterations = 3
                )
            stdout_psi_uni, stderr_psi_uni = psiblast_uniprot_cline()
            psiblast_pdb_cline = NcbipsiblastCommandline(
                db = 'db/PDBdb', out = "psiblast_pdb_results.xml", outfmt = 5,
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
                if len(hits_dict) == 7:
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
    """
    Downloads a list of PDB files and stores them in the structures directory

    Input: List of PDB codes of the proteins

    Return: Directory called "structures" with the PDB files in it in the format: "structures/pdb{code}.ent"
    """
    from Bio.PDB.PDBList import PDBList
    ## Get PDB files
    r = PDBList()  # Object instance
    r.download_pdb_files(pdb_codes = pdb_codes_list , file_format = "pdb",
                            pdir = "structures/", overwrite=True)   # creates the directory if do not exists

def pdb_quality_filter(protein_codes, hits_dict, resolution = 2):
    # Needs the protein_codes and hits_dict from the get_pdb_homologs function
    """
    Takes the files from the "structures/" directory, parses them to obtain the resolution and returns the three
    hits with the highest e-value and a resolution higher than 2 Angstroms.

    Input: List of protein codes from the BLAST and dictionary of the chains for each hit. Optional: Resolution

    Return: List of three protein PDB codes and their respective chains
    """
    from Bio.PDB import PDBParser
    import pathlib

    res_dict = {}
    for pdb_file in pathlib.Path("structures/").iterdir():
        # Read pdb file
        with open(pdb_file, "r") as fh:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure(fh.name[14:18], fh)
            res_dict[fh.name[14:18]] = structure.header['resolution']
    # Order the dictionary according to the highest e-value order
    ordered_res_dict = {k: res_dict[k] for k in protein_codes}

    # Now, we will keep just the three hits with the highest e-value and a resolution lower than 2
    final_proteins = []
    final_chains = []

    for key in ordered_res_dict:
        if len(final_proteins) == 3:
            break
        elif ordered_res_dict[key] <= resolution:
            final_proteins.append(key)

    for code in final_proteins:
        final_chains.append(hits_dict[code])

    return final_proteins, final_chains

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
        - List containing the residue length associated to the pbd sequence used in the alignment
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
    """Obtain the normalized b-factors of c-alpha atoms
    FORMULA!!!

    Return: list of normalized b-factors
    """
    import statistics
    n_bfactor = list(map(lambda bf:
                (bf-statistics.mean(bfactor))/statistics.stdev(bfactor), bfactor))
    return n_bfactor

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

    ## Get the b-factor only for regions of similarity
    all_sim_bfactors = list()

    for val in list(zip(all_norm_bfactors, all_sim)):
        norm_bfactors = val[0]
        similarities  = val[1]

        sim_bfactors = [norm_bfactors[i] for i in similarities]
        all_sim_bfactors.append(sim_bfactors)

    ## Compute the mean
    final_bfactors = [sum(sim_bfactors)/len(all_sim_bfactors) for sim_bfactors in zip(*all_sim_bfactors)]
