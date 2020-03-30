#!/bin/bash

module load ANACONDA/2018.12_py3
source activate biobb

#COMPSs environment
export COMPSS_PYTHON_VERSION=none
#COMPSs release
module load COMPSs/2.6

#Singularity
module load singularity

#GROMACS 2019
module load intel/2018.4 impi/2018.4 mkl/2018.4 gromacs/2019.1


enqueue_compss --num_nodes=16 --exec_time=120 --network=ethernet --worker_working_dir=$PWD --master_working_dir=$PWD --base_log_dir=$PWD \
 --qos=debug  --job_name=T790M_IRE --worker_in_master_cpus=48 --jvm_workers_opts="-Dcompss.worker.removeWD=false" \
 pmxLig.py --config pmxLig.yaml

