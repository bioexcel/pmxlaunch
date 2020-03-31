#!/bin/bash

module purge

module load ANACONDA/2018.12_py3
source activate biobb

# COMPSs environment
export COMPSS_PYTHON_VERSION=none
# COMPSs release
module load COMPSs/2.6.1

# Singularity
module load singularity

#GROMACS 2019
module load intel/2018.4 impi/2018.4 mkl/2018.4 gromacs/2019.1

#Permissions for everyone
umask ugo+rwx

enqueue_compss -d --num_nodes=3 --exec_time=120 \
	--base_log_dir=$PWD \
	--network=ethernet --qos=debug \
	--pythonpath=$PWD/wf_defs $PWD/wf_defs/main.py 
