#/bin/bash

export PSRHOME=/fred/oz002/psrhome
source ${PSRHOME}/scripts/psrhome.sh

# Unload psrhome Python modules
module unload anaconda2/5.1.0
module load python/3.6.4
module load numpy/1.16.3-python-3.6.4
module load scipy/1.3.0-python-3.6.4
module load psrchive/b7ad99c37-python-3.6.4
module load astropy/3.1.2-python-3.6.4
