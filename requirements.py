import subprocess
import os

# Blastp and psiblast: https://www.ncbi.nlm.nih.gov/books/NBK569861/
subprocess.run(["wget", "https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.13.0+-x64-linux.tar.gz"], check=True)
subprocess.run(["tar", "zxvpf", "ncbi-blast-2.13.0+-x64-linux.tar.gz"], check=True)
subprocess.run(["export", "PATH=$PWD/ncbi-blast-2.13.0+/bin:$PATH"], check=True) # temporary solution. profile must be modified

# Tcoffe: https://www.tcoffee.org/Projects/tcoffee/workshops/tcoffeetutorials/installation.html
subprocess.run(["sudo", "apt-get", "install", "t-coffee"], check=True) # sudo apt-get install t-coffee

# Dssp: https://github.com/cmbi/dssp
subprocess.run(["sudo", "apt-get", "install", "libboost-all-dev"], check=True)
subprocess.run(["wget", "https://github.com/cmbi/dssp/archive/refs/tags/2.3.0.tar.gz"], check=True)
subprocess.run(["tar", "-zxvf", "2.3.0.tar.gz"], check=True)
os.chdir("dssp-2.3.0/")

subprocess.run(["./autogen.sh"], check=True)
subprocess.run(["./configure"], check=True)
subprocess.run(["make"], check=True)
subprocess.run(["make", "mkdssp"], check=True)
subprocess.run(["sudo","make", "install"], check=True)
