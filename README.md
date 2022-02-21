# Setup

Please, first of all, download this repository locally:
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
