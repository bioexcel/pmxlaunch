#!/usr/bin/env python3

import os
import sys
import zipfile
import time
import argparse

# biobb common modules
from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu

# pycompss: biobb analysis modules
from biobb_adapters.pycompss.biobb_analysis.gromacs.gmx_image_pc import gmx_image_pc

def main(config, system=None):
    conf = settings.ConfReader(config, system)
    global_log, _ = fu.get_logs(path=conf.get_working_dir_path(), light_format=True)

    for mutation in conf.properties['mutations']:
        for ligand in conf.properties['ligands']:

            traj_code = "_".join([mutation, ligand])
            traj_code_path = os.path.join(mutation, ligand)
            ensemble_prop = conf.get_prop_dic(prefix=traj_code_path)
            ensemble_paths = conf.get_paths_dic(prefix=traj_code_path)

            # Image
            global_log.info("Imaging " + traj_code + " trajectory to remove PBC issues")
            # /gpfs/projects/bsc23/bsc23513/BioExcel/BioExcel_EGFR_pmx/MDs/WT/WT_apo_md_FULL.xtc
            folder = conf.properties['workdir']+"/"+mutation
            traj = folder + "/" +  traj_code + "_md_FULL.xtc"
            top = folder + "/" + traj_code + "_md_1.tpr"
            out = folder + "/" +  traj_code + "_md_FULL.imaged.xtc"
            global_log.info("Traj in:" + traj)
            global_log.info("Top in:" + top)
            global_log.info("Traj out:" + out)
            ensemble_paths['image']['input_traj_path'] = traj
            ensemble_paths['image']['input_top_path'] = top
            ensemble_paths['image']['output_traj_path'] = out
            gmx_image_pc(**ensemble_paths["image"], properties=ensemble_prop["image"])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Biobb trajectory imaging with GROMACS")
    parser.add_argument('--config', required=True)
    parser.add_argument('--system', required=False)
    args = parser.parse_args()
    main(args.config, args.system)
