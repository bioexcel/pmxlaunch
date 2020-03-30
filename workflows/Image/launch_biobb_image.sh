#!/bin/bash

source env2.sh

enqueue_compss --num_nodes=4 --exec_time=120 --network=ethernet --worker_working_dir=$PWD --master_working_dir=$PWD --base_log_dir=$PWD \
 --qos=debug -d -g --worker_in_master_cpus=48 --jvm_workers_opts="-Dcompss.worker.removeWD=false" \
 /gpfs/projects/bsc23/bsc23513/BioExcel/BioExcel_EGFR_pmx/MDs/Adam/Wfs/Image/biobb_image.py --config /gpfs/projects/bsc23/bsc23513/BioExcel/BioExcel_EGFR_pmx/MDs/Adam/Wfs/Image/biobb_image.reduced.yaml

