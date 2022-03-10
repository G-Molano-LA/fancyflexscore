# Setup

Please, first of all, download this repository locally :
```
git clone https://github.com/G-Molano-LA/structured-python.git
```

# Project Objective

Programming a standalone program for solving a specific problems.


 **Objective**: To develop a flexibility score for proteins
 **Input**: protein sequence or protein family
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

## Baldo's remarks
You can use the b-factor to calculate the protein flexibility, however, this parameter is not totally correct as it contains cristall vibrations.
