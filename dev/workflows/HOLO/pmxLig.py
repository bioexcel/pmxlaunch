#!/usr/bin/env python3

import os
import zipfile
import sys
import time
import argparse
from pycompss.api.api import compss_wait_on_file

# biobb common modules
from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu

# pycompss: biobb pmx modules
from biobb_adapters.pycompss.biobb_pmx.pmx.mutate_pc import mutate_pc
from biobb_adapters.pycompss.biobb_pmx.pmx.gentop_pc import gentop_pc
from biobb_adapters.pycompss.biobb_pmx.pmx.analyse_pc import analyse_pc

# pycompss: biobb md modules
from biobb_adapters.pycompss.biobb_md.gromacs.pdb2gmx_pc import pdb2gmx_pc
from biobb_adapters.pycompss.biobb_md.gromacs.make_ndx_pc import make_ndx_pc
from biobb_adapters.pycompss.biobb_md.gromacs.grompp_pc import grompp_pc
from biobb_adapters.pycompss.biobb_md.gromacs.grompp_ndx_pc import grompp_ndx_pc
from biobb_adapters.pycompss.biobb_md.gromacs.mdrun_pc import mdrun_pc
from biobb_adapters.pycompss.biobb_md.gromacs.mdrun_dhdl_pc import mdrun_dhdl_pc
from biobb_adapters.pycompss.biobb_md.gromacs_extra.append_ligand_pc import append_ligand_pc

# pycompss: biobb analysis modules
from biobb_adapters.pycompss.biobb_analysis.gromacs.gmx_image_pc import gmx_image_pc
from biobb_adapters.pycompss.biobb_analysis.gromacs.gmx_trjconv_str_ens_pc import gmx_trjconv_str_ens_pc

# pycompss: biobb structure utils modules
from biobb_adapters.pycompss.biobb_structure_utils.utils.extract_atoms_pc import extract_atoms_pc
from biobb_adapters.pycompss.biobb_structure_utils.utils.remove_ligand_pc import remove_ligand_pc
from biobb_adapters.pycompss.biobb_structure_utils.utils.sort_gro_residues_pc import sort_gro_residues_pc


def main(config, imaged_traj_available):
    start_time = time.time()
    conf = settings.ConfReader(config)
    global_log, _ = fu.get_logs(path=conf.get_working_dir_path(), light_format=True)
    global_prop = conf.get_prop_dic()
    global_paths = conf.get_paths_dic()

    dhdl_paths_listA = []
    dhdl_paths_listB = []
    for ensemble, mutation in conf.properties['mutations'].items():

        ensemble_prop = conf.get_prop_dic(prefix=ensemble)
        ensemble_paths = conf.get_paths_dic(prefix=ensemble)

        if not imaged_traj_available:
            # step0_image There isn't an imaged trajectory
            global_log.info(ensemble + " Step 0: gmx image: Imaging trajectories to remove PBC issues")
            ensemble_paths['step0_image']['input_top_path'] = conf.properties['input_trajs'][ensemble]['input_tpr_path']
            ensemble_paths['step0_image']['input_traj_path'] = conf.properties['input_trajs'][ensemble]['input_traj_path']
            gmx_image_pc(**ensemble_paths["step0_image"], properties=ensemble_prop["step0_image"])

        else:
            # An imaged trajectory is available
            ensemble_paths['step1_trjconv_'+ensemble]['input_traj_path'] = conf.properties['input_trajs'][ensemble]['input_traj_path']

        # step1_trjconv
        global_log.info(ensemble + " Step 1: gmx trjconv: Extract snapshots from equilibrium trajectories")
        ensemble_paths['step1_trjconv_'+ensemble]['input_top_path'] = conf.properties['input_trajs'][ensemble]['input_tpr_path']
        gmx_trjconv_str_ens_pc(**ensemble_paths['step1_trjconv_'+ensemble], properties=ensemble_prop['step1_trjconv_'+ensemble])

    for ensemble, mutation in conf.properties['mutations'].items():
        ensemble_prop = conf.get_prop_dic(prefix=ensemble)
        ensemble_paths = conf.get_paths_dic(prefix=ensemble)
        compss_wait_on_file(ensemble_paths['step1_trjconv_'+ensemble]["output_str_ens_path"])

        with zipfile.ZipFile(ensemble_paths['step1_trjconv_'+ensemble]["output_str_ens_path"]) as zip_f:
            unique_dir = os.path.abspath(fu.create_unique_dir(prefix=ensemble_prop['step1_trjconv_'+ensemble]['working_dir_path']+'/'+ensemble+'/'))
            zip_f.extractall(unique_dir)
            state_pdb_list = [os.path.join(unique_dir, name)for name in zip_f.namelist()]

        for pdb_path in state_pdb_list:
            pdb_name = os.path.splitext(os.path.basename(pdb_path))[0]
            prop = conf.get_prop_dic(prefix=os.path.join(ensemble, pdb_name))
            paths = conf.get_paths_dic(prefix=os.path.join(ensemble, pdb_name))

            # step2_pmx_mutate
            global_log.info(ensemble+" "+pdb_name+" Step 2: pmx mutate: Generate Hybrid Structure")
            paths['step2_pmx_mutate']['input_structure_path'] = pdb_path
            prop['step2_pmx_mutate']['mutation_list'] = mutation
            mutate_pc(**paths["step2_pmx_mutate"], properties=prop["step2_pmx_mutate"])

            # step3_check_dummies
            global_log.info(ensemble + " " + pdb_name + " Step 3 Check for dummy atoms")
            extract_atoms_pc(**paths['step3_check_dummies'], properties=prop['step3_check_dummies'])

            compss_wait_on_file(paths['step3_check_dummies']['output_structure_path'])
            try:
                dummy = bool(os.path.getsize(paths['step3_check_dummies']['output_structure_path']))
            except:
                global_log.info("Error ocurred while checking the file containing dummy atoms")
                global_log.info(sys.exc_info()[0])

            # step4_remove_ligand
            global_log.info(ensemble + " " + pdb_name + " Step 4 Remove ligand")
            remove_ligand_pc(**paths['step4_remove_ligand'], properties=prop['step4_remove_ligand'])

            # step5_gmx_pdb2gmx
            global_log.info(ensemble + " " + pdb_name + " Step 5: gmx pdb2gmx: Generate Topology")
            pdb2gmx_pc(**paths["step5_gmx_pdb2gmx"], properties=prop["step5_gmx_pdb2gmx"])

            # step6_sort_gro
            global_log.info(ensemble + " " + pdb_name + " Step 6 Sort gro residues")
            sort_gro_residues_pc(**paths['step6_sort_gro'], properties=prop['step6_sort_gro'])

            # step7_lig_gmx_appendLigand
            global_log.info(ensemble + " " + pdb_name + " Step 7_lig: gmx appendLigand: Append a ligand to a GROMACS topology")
            append_ligand_pc(**paths["step7_lig_gmx_appendLigand"], properties=prop["step7_lig_gmx_appendLigand"])

            # step8_pmx_gentop
            global_log.info(ensemble + " " + pdb_name + " Step 8: pmx gentop: Generate Hybrid Topology")
            gentop_pc(**paths["step8_pmx_gentop"], properties=prop["step8_pmx_gentop"])

            if not dummy:
                paths['step12_gmx_grompp']['input_gro_path'] = paths['step6_sort_gro']['output_gro_path']
            else:
                # step9_gmx_makendx
                global_log.info(ensemble + " " + pdb_name + " Step 9 (Dummies): gmx make_ndx: Generate Gromacs Index file to select atoms to freeze")
                make_ndx_pc(**paths["step9_gmx_makendx"], properties=prop["step9_gmx_makendx"])

                # step10_gmx_grompp
                global_log.info(ensemble + " " + pdb_name + " Step 10 (Dummies): gmx grompp: Creating portable binary run file for energy minimization")
                grompp_ndx_pc(**paths["step10_gmx_grompp"], properties=prop["step10_gmx_grompp"])

                # step11_gmx_mdrun
                global_log.info(ensemble + " " + pdb_name + " Step 11 (Dummies): gmx mdrun: Running energy minimization")
                mdrun_pc(**paths["step11_gmx_mdrun"], properties=prop["step11_gmx_mdrun"])

            # step12_gmx_grompp
            global_log.info(ensemble + " " + pdb_name + " Step 12: gmx grompp: Creating portable binary run file for system equilibration")
            grompp_pc(**paths["step12_gmx_grompp"], properties=prop["step12_gmx_grompp"])

            # step13_gmx_mdrun
            global_log.info(ensemble + " " + pdb_name + " Step 13: gmx mdrun: Running system equilibration")
            mdrun_pc(**paths["step13_gmx_mdrun"], properties=prop["step13_gmx_mdrun"])

            # step14_gmx_grompp
            global_log.info(ensemble + " " + pdb_name + " Step 14: Creating portable binary run file for thermodynamic integration (ti)")
            grompp_pc(**paths["step14_gmx_grompp"], properties=prop["step14_gmx_grompp"])

            # step15_gmx_mdrun
            global_log.info(ensemble + " " + pdb_name + " Step 15: gmx mdrun: Running thermodynamic integration")
            mdrun_dhdl_pc(**paths["step15_gmx_mdrun"], properties=prop["step15_gmx_mdrun"])

            if ensemble == "stateA":
                dhdl_paths_listA.append(paths["step15_gmx_mdrun"]["output_dhdl_path"])
            elif ensemble == "stateB":
                dhdl_paths_listB.append(paths["step15_gmx_mdrun"]["output_dhdl_path"])

    # Creating zip file containing all the dhdl files
    dhdlA_path = os.path.join(global_prop["step16_pmx_analyse"]["working_dir_path"],'dhdlA.zip')
    dhdlB_path = os.path.join(global_prop["step16_pmx_analyse"]["working_dir_path"],'dhdlB.zip')
    for dhdl_file in dhdl_paths_listA:
        compss_wait_on_file(dhdl_file)
    for dhdl_file in dhdl_paths_listB:
        compss_wait_on_file(dhdl_file)
    fu.zip_list(dhdlA_path, dhdl_paths_listA)
    fu.zip_list(dhdlB_path, dhdl_paths_listB)

    # step16_pmx_analyse
    global_log.info(ensemble+" Step 16: pmx analyse: Calculate free energies from fast growth thermodynamic integration simulations")
    global_paths["step16_pmx_analyse"]["input_a_xvg_zip_path"]=dhdlA_path
    global_paths["step16_pmx_analyse"]["input_b_xvg_zip_path"]=dhdlB_path
    analyse_pc(**global_paths["step16_pmx_analyse"], properties=global_prop["step16_pmx_analyse"])

    elapsed_time = time.time() - start_time
    global_log.info('')
    global_log.info('')
    global_log.info('Execution successful: ')
    global_log.info('  Workflow_path: %s' % conf.get_working_dir_path())
    global_log.info('  Config File: %s' % config)
    global_log.info('')
    global_log.info('Elapsed time: %.1f minutes' % (elapsed_time/60))
    global_log.info('')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Based on the official PMX tutorial")
    parser.add_argument('--config', required=True)
    parser.add_argument('--imaged_traj_available', required=False, action='store_true')
    args = parser.parse_args()
    main(config=args.config, imaged_traj_available=args.imaged_traj_available)
