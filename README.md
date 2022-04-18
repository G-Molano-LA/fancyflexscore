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

# TO DO LIST
8. Do documentation:
  -  [ ] Configure logging for all the modules
  -  [ ] Tutorial with examples on how to use the program
  -  [ ] Analysis of examples of 4 cases. The results analyis should contain at least two cases of the following list, and two other cases (not necessarily in this list).


# LIMITATIONS
- Our approach do not take into account unkwnon regions as we do not obtain b-factors. To supply this missing values we thought to put the mean of b-factor, but, as we standarize this value later, we decided that makes no sense to put the mean (as gives 0 in standarization).

# Project Objective
Programming a standalone program for solving a specific problems.

 **Objective**: To develop a flexibility score for proteins  
 **Input**: protein sequence or protein family. The input should be a sequence in FASTA format, not a uniprot code.
 **Output**:  
 - Flexibility score for each aminoacid in the protein sequence
  - Parseable text output files
  - Graphical representation of the scores


PYT evaluation will take into account:
- Program structure (tasks, classes, modules…)
- Reusability, use of libraries.
- Others

- Prepare the structure of the program:
  - Prepare the package directory structure
  - Prepare the setup.py file

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

## Articles of interest

[Functional aspects of protein flexibility](https://link.springer.com/article/10.1007/s00018-009-0014-6)  
[Utility of B-Factors in Protein Science](https://pubs.acs.org/doi/10.1021/acs.chemrev.8b00290)  
[Flexibility is a mechanical determinant of antimicrobial activity for amphipathic cationic α-helical antimicrobial peptides](https://www.sciencedirect.com/science/article/pii/S000527361300206X?via%3Dihub#s0010)  
[Flexibility Function - source code](https://home.cc.umanitoba.ca/~psgendb/doc/local/biopython-1.64.old/Bio/SeqUtils/ProtParam.py)
