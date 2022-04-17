
## Uniprot codes
import Bio.SwissProt as sp
## To download files from web
import requests
## Catch the stdoute
from io import StringIO
import sys
## Regex
import re
import logging

# General configuration of logging module. More info: https://docs.python.org/3/howto/logging.html#logging-basic-tutorial
logging.basicConfig(filename='flexscore.log', encoding = 'utf-8', level=20) # level logging.info


class capturing_output(list):
    """Captures STDOUT by creating a new STDOUT

    The special methods (__enter__, __exit__) allows create class instances by
    using the 'with' statement.

    Example:
        with capturing_output() as out:
            do_something(my_object)
        out.getvalue()
    """

    def __enter__(self):
        # Save the current stdout in a class atttribtue
        self._stdout = sys.stdout
        # and create a new one, save it as a class attribute and replace the previous one
        sys.stdout = self._stringio = StringIO()
        return self
    def __exit__(self, *args):
        # return a list object
        self.extend(self._stringio.getvalue().splitlines())
        # Delete created STDOUT to free up some memory
        del self._stringio
        # Restore the initial stdout
        sys.stdout = self._stdout

# Obtain PDB code for target protein, if any
## Obtain uniprot file
url       = f"https://www.uniprot.org/uniprot/{target}.txt"
localfile = f"Uniprot/{target}.txt"
get_stream_data(url, localfile)
## Read uniprot entrance
with open(localfile) as fh:
    uniprot_record = sp.read(fh)
    print(uniprot_record.references)
    print(uniprot_record.cross_references)

# If there is not pdb file, get from alpha fold
if "doesn't exists" in captured_out[1]:
    """
    The download_pdb_files function always tries to download the pdb file,
    printing:
        "Downloading PDB structure 'pdb_code'...""
    We the structure is not found, prints:
        "Desired structure doesn't exists"
    Therefore, if there is no structure for the target, we will download the structure
    from alphafold
    """
    # Get the protein_code
    logging.info(f"{target} has not PDB file. Getting structure from AlphaFold database...")
    # Obtain pdb file from AlphaFold
    ## Define the remote file to retrieve
    url = f"https://alphafold.ebi.ac.uk/files/AF-{target}-F1-model_v2.pdb"
    # Define the local filename to save data
    local_file = f'structures/{target}.pdb'
    # Make http request for remote file data
    get_stream_data(url, localfile)

def get_stream_data(url, outfile_path):
    """Download data from websites
    """
    import requests

    # Make http request for remote file data
    target_data = requests.get(url)
    # Save file data to local copy
    with open(outfile_path, 'wb')as file:
        file.write(target_data.content)

from pathlib import Path
from glob import glob
# Change the extension of pdb '.ent' to '.pdb'
files = glob("structures/*.ent")
for name in files:
    p = Path(name)
    s = p.rename(p.with_suffix('.pdb'))

pdb_outfiles = [f"structures/pdb{code}.pdb" for code in final_protein_codes]

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
    all_rest = list()
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
