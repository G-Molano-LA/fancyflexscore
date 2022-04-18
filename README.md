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
For running the program, the user should execute:
```
python3 fancyflexprot -i input_file -o output_file_name -ws window_size -v
```
