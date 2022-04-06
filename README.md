# Setup

Please, first of all, download this repository locally :
```
git clone https://github.com/G-Molano-LA/structured-python.git
```
# KEY IDEA

The protein flexibility depends on several parameters. There are some parameter that are indicators of that:
- B-factor
- Secondary structure  
- Hidrofobicity

Then, the key idea is to use this 3 parameters to calculate our flexibility score.

To begin with, we have a protein sequence as a input that has no known pdb file. So.
1. Search candidates that has similar structure with our target protein -> Pairwise aln (BLAST, jackhmmer or similar approaches)
2. With these candidates, obtain regions of similarity -> Multiple **STRUCTURAL** aln
3. Obtain b-factors of regions of similarity (perform an average of the candidates).
4. Compute b-factor calculation:
  1. Obtain the b-factor associated to the alfa-carbon of each aminoacid.
  2. Standarize b-factors.
5. Compute a flexibility score for **each aminoacid** based on:
  - B-factors
  - Secondary structure restraints (Can be calculated with the Bio.PDB.DSSP module)
  - Hidrofobicity (Can be calculated with the Bio.PDB.DSSP module)
6. Define a threshold to decide if the score for each aminoacid is flexible or not (1, 0).
7. Scan the sequence by triplets to check for different state residues in the middle of the triplet to correct it.
8. Sum all 0,1s and normalize by the number of aminoacids (average).
9. As a final output we can give:
  - Total protein score
  - Score associated to each aminoacid

# TO DO LIST

1. Search candidates that has similar structure with our target protein: GERARD
  INPUT: Fasta Sequence of the target protein
  OUTPUT: PDB codes
  - [ ] Decide which workflow is better to follow to perform a pairwise aln (blast, psiblast, jackhammer...)
  - [ ] Implement the pairwise aln in python
2. Obtain the PDB model of the target protein: ALE
  INPUT: PDB codes
  OUTPUT: PDB files
  - [ ] Obtain PDB from alpha-fold is the sequence is not known
  - [ ] Obtain PDB from database if the protein is known
3. With these candidates, obtain regions of similarity: IRIA
  - [ ] Decide which multiple **STRUCTURAL** aln (MSA) approach to use
  - [ ] Implement the MSA in python
4. Obtain b-factors of regions of similarity (implement in python):
  - [ ] Identify regions of similarity in the MSA
  - [ ] Take the b-factor of this regions of similarity:
    - If there there is only one candidate sharing this region: take the alfa-carbon b-factor of the aminoacids
    - If there is more than one candidate sharing this region: take the mean of the alfa-carbon b-factor of the aminoacid
  - [ ] Standaritze values
5. Compute values based on:
  - [ ] Structural restraints.
  - [ ] Hidrofobicity.
6. Give weights to the different parameters:
    - [ ] B-factors ex. 0.33
    - [ ] Secondary structure restraints ex. 0.33
    - [ ] Hidrofobicity ex. 0.33  
7. Obtain a final flexibility score:
  - [ ] Weighted score for each aminoacid
  - [ ] Define a threshold to decide if the weighted score is flexible or not = obtain 0 and 1s.
  - [ ] Obtain the final score for each aminoacid by applying the threshold. (airi formula)
    $$ f(x) = $$
  - [ ] Sum all 0,1s and normalize by the number of aminoacids (average). (airi formula)
    $$ \sigma = $$
8. Do documentation:
  -  [ ] Configure logging for all the modules
  -  [ ] Tutorial with examples on how to use the program
  -  [ ] Analysis of examples of 4 cases. The results analyis should contain at least two cases of the following list, and two other cases (not necessarily in this list).


# LIMITATIONS

- Our approach do not take into account unkwnon regions as we do not obtain b-factors. To supply this missing values we thought to put the mean of b-factor, but, as we standarize this value later, we decided that makes no sense to put the mean (as gives 0 in standarization).

# DEPENDENCIES
- ICM installed

# REMARKS  
- Baldo's said : "You can use the b-factor to calculate the protein flexibility, however, this parameter is not totally correct as it contains cristall vibrations."
- Gery said: b-factor del alpha fold = bad predictions (all have high values)  



# Project Objective

Programming a standalone program for solving a specific problems.


 **Objective**: To develop a flexibility score for proteins  
 **Input**: protein sequence or protein family. The input should be a sequence in FASTA format, not a uniprot code.
 **Output**:  
 - Flexibility score for each aminoacid in the protein sequence
  - Parseable text output files
  - Graphical representation of the scores

Proposed steps:
- Multiple Sequence Alignment
- Structural Alignment
- Develop a flexibility score analyzing the alignments
- Graphical  representation of the flexibility scores
- Results analyzed for 2 specified families
- Results analyzed for 2 other families

PYT evaluation will take into account:
- Program structure (tasks, classes, modules…)
- Reusability, use of libraries.
- Others

# PYT Session 11 Objective
Prepare the input and output interface for the SBI-PYT project.
- When your project is executed as a standalone program, it will be
expected to have some input and output arguments.

# PYT Session 12 Objective
- Install third-party packages
  - Biopython, scipy, numpy, pandas, matplotlib, seaborn
- Prepare the structure of the program:
  - Prepare the package directory structure
  - Prepare the setup.py file

## The teacher has recommended to use a virtual environment to test the dependencies of the program. To use a virtual environment, se can use:
```
sudo apt-get install python3.9-venv
pip install pyvenv
python3 -m venv structural_py
```
Then, you will have a directory named `structured_py` which contains all the programs of the virtual env in the `bin` directory. Then, if you have to install something inside this virtual env, use the pip3 programm located at the `bin` directory. The installed programs will be saved on the `lib` directory.

Note: I have included `structured_py/` in the `.gitignore` file, in order to avoid the upload of the virtual env in the repository.

## Writting the setup Script
```
from distutils.core import setup

setup(name='Structured Python',
  version='1.0',
  description='Flexibility score for proteins',
  author='Iria Pose, Gerard Romero Sola, Leidy Alejandra Gonzalez Molano',
  author_email='javigarcia@python.org',
  url='',
  packages=['', ''],
)
```
Thanks!

### Possible approach

The **B-factor is** a cristallography parameter to determine flexibility. It is based on the movement of the different residues when applying temperature (it's not exactly like that, but you get the concept). So, it could be a good parameter to get the score, but some more restrains could be applied in order to make the score even better: **type of secondary structure, hydrophobicity regions and aminoacid preference** (how aminoacids interact with the other aminoacids in the region). However, as a first approach it could be enough.
So, this parameter is obtained from crystallography approaches and it's found in the last numeric column of the PDB file By knowing the ranges of this parameter (15-30 as rigid and above for flexible regions) we can make some score associations. Taking this into account, we have an open source app which is AlphaFold2, which provides a PDB file for the predicted protein where the B-factor appears.
I realized the prediction of the B-factor is quite bad, since it gives very big values, so maybe we should scale them or do some kind of cross-validation with the structure prediction found in the mmCIF file. Or just use this structure prediction to assess scoring.
In case of finally using the b-factor somehow, we can normalize it by z-scores which is done in the following reference:
https://www.blopig.com/blog/2015/08/using-b-factors-to-assess-flexibility/

We also have servers or applications for flexibility prediction like MEDUSA: https://www.dsimb.inserm.fr/MEDUSA/index.html, which can be useful for assessing if our scoring is properly done.

As an example of other year's project:
https://github.com/martaferri/SBI_Project

#### Python's teacher remarks

Since implementing Alphafold in python is computationally very expensive, we can extract pdb files from the AlphaFold's already predicted repository or from the PDB if the structure is known.
