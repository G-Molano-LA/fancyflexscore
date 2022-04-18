# The FancyFlexScore

**Here, you will find a short review about FancyFlexScore. Further information can be found in report pdf.**  
The main purpose of this project is to provide an assessment of the flexibility of a protein given its sequence. In detail, the function outputs a text file containing the flexibility scores per amino acid in addition to a graphical representation of these scores along the sequence. Furthermore, both outputs include complementary information regarding the hydrophobicity and secondary structures of the protein, to provide a comparison with parameters theoretically related to flexibility.

# Prerequisites
Python Packages:
- biopython >= v.1.79
- matplotlib >= v.3.3.4
- statistics >=
- requests >= v.2.25.1
- pandas >=
- seaborn >= v.0.11.2

The above packages can be easily installed by executing the `requirements.txt` file:
```
pip install -r requirements.txt
```
External Programs. The following programs must be on the `$PATH`:
- [NCBI programs](https://www.ncbi.nlm.nih.gov/books/NBK569861/) (Blastp, Psiblast and makeblastdb) >= v.2.10.1+:  Pairwise alignment and database creation.
- [T-Coffee](https://www.tcoffee.org/Projects/tcoffee/workshops/tcoffeetutorials/installation.html) >= v.13.45.0: Multiple sequence alignment
- [DSSP](https://github.com/cmbi/dssp) >= v.2.3.0 (mkdssp): Secondary Structure assignment.

The installation of the above programs can be automatically handle by executing the `requirements.py` file (linux users):
```
 python3 requirements.py
```

# Running the program
The user can obtain information about the program and the arguments that should write by:
```
python3 fancyflexprot -h
```

For running the program, the user should execute:
```
python3 fancyflexprot -i input_file -o output_file_name -ws window_size -v
```
The arguments are:
- -i INFILE: Input FASTA file with one unique record containing the identifier and sequence of the protein of interest. The identifier must be a uniprot ID. Chain must be specified by: >uniprotID:chain. This argument is required.  
- -o OUTFILE: Output filename without extensions. It generates two files: filename\_results.csv and filename\_visualization.pdf (default: results). Two more files are also generated: a PDB for homologs (located at structures/directory) and a MSA aln of target protein and PDB homologs ((located at structures/directory)). THis argument is optional.  
- -v verbose: To have a follow up of the program (default: False). This argument is optional.    
- -ws window\_size: To change the window size when computing the flexibility score. The range of accepted values: {1,3,5,7,9,11,13}. This argument is optional (default: 7).  
