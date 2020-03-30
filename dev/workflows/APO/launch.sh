#!/bin/bash

source env.sh

enqueue_compss --num_nodes=16 --exec_time=120 --network=ethernet --worker_working_dir=$PWD --master_working_dir=$PWD --base_log_dir=$PWD \
 --qos=debug --job_name=T790M_apo --worker_in_master_cpus=48 --jvm_workers_opts="-Dcompss.worker.removeWD=false" \
 pmx_apo.py --config pmx_apo.yaml

