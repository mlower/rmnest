#/bin/bash

export PSRHOME=/fred/oz002/psrhome
source ${PSRHOME}/scripts/psrhome.sh

# Unload psrhome Python modules
module unload anaconda2/5.1.0
module load python/3.7.4
module load numpy/1.18.2-python-3.7.4
module laod scipy/1.6.0-python-3.7.4
module load psrchive/c216582a0-gcc-9.2.0-python-3.7.4
module load astropy/4.0.1-python-3.7.4
