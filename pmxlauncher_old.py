#!/usr/bin/env python3

from pathlib import Path
import argparse
from Bio.PDB.Polypeptide import three_to_one
from Bio.PDB.Polypeptide import one_to_three
import re
import shutil
import subprocess
import oyaml as yaml

RESIDUE_NUMBER_OFFSET = 693


def create_biobb_pth_file(file_path):
    with open(file_path, 'w') as pth_file:
        pth_file.write("/home/bsc23/bsc23210/macshare/biobb_common \n")
        pth_file.write("/home/bsc23/bsc23210/macshare/biobb_md \n")
        pth_file.write("/home/bsc23/bsc23210/macshare/biobb_pmx \n")
        pth_file.write("/home/bsc23/bsc23210/macshare/biobb_analysis \n")
        pth_file.write("/home/bsc23/bsc23210/macshare/biobb_adapters \n")
        pth_file.write("/home/bsc23/bsc23210/macshare/biobb_wf_pmxligand \n")
        pth_file.write("/home/bsc23/bsc23210/macshare/biobb_structure_utils")


#def get_mutation_dict(mutation):
#    pattern = re.compile(r"(?P<wt>[a-zA-Z]{3})(?P<resnum>\d+)(?P<mt>[a-zA-Z]{3})")
#    return pattern.match(mutation.strip()).groupdict()

def get_mutation_dict(mutation):
    if mutation.strip()[1].isdigit():
        pattern = re.compile(r"(?P<wt>[a-zA-Z]{1})(?P<resnum>\d+)(?P<mt>[a-zA-Z]{1})")
        mut_dict = pattern.match(mutation.strip()).groupdict()
        mut_dict['wt'] = one_to_three(mut_dict['wt'].upper())
        mut_dict['mt'] = one_to_three(mut_dict['mt'].upper())
    else:
        pattern = re.compile(r"(?P<wt>[a-zA-Z]{3})(?P<resnum>\d+)(?P<mt>[a-zA-Z]{3})")
        mut_dict = pattern.match(mutation.strip()).groupdict()
    return mut_dict


def three_to_one_mutation(mutation):
    mut_dict = get_mutation_dict(mutation)
    return f"{three_to_one(mut_dict.get('wt').upper())}{mut_dict.get('resnum')}{three_to_one(mut_dict.get('mt').upper())}"


def get_input_tpr_xtc_paths(name, replica, ligand, base_dir):
    modifier = ''
    if '_' in ligand:
        ligand, modifier = tuple(ligand.split('_'))

    str_replica = 'x' + str(replica)
    if replica == 1:
        str_replica = ''
    if name.lower() == 'wt':
        traj_dir_name = 'WT' + str_replica
    else:
        traj_dir_name = f"{three_to_one_mutation(name)}{str_replica}"

    traj_ligand_name = f'{traj_dir_name}_{ligand.lower()}'
    if modifier:
        traj_ligand_name = f'{traj_ligand_name}_{modifier.upper()}'

    traj_dir_path = base_dir.joinpath('MDs', traj_dir_name)
    traj_tpr_path = traj_dir_path.joinpath(f'{traj_ligand_name}_md_1.tpr')
    traj_xtc_path = traj_dir_path.joinpath(f'{traj_ligand_name}_md_FULL.imaged.xtc')
    if not traj_xtc_path.exists():
        traj_xtc_path = traj_dir_path.joinpath(f'{traj_ligand_name}_md_FULL.xtc')
    return traj_tpr_path, traj_xtc_path


def is_imaged_traj_available(wt_traj, mut_traj, pattern='imaged'):
    return pattern in str(wt_traj) and pattern in str(mut_traj)


def get_template_config_dict(config_yaml_path):
    with open(config_yaml_path) as config_yaml_file:
        return yaml.safe_load(config_yaml_file)


def launch(mutation, pmx_resnum, ligand, wt_replica, mut_replica, queue, num_nodes, compss_version, fe_nsteps, trjconv_skip, base_dir, compss_debug, time):
    if pmx_resnum == 0:
        pmx_resnum = int(get_mutation_dict(mutation)['resnum']) - RESIDUE_NUMBER_OFFSET

    base_dir = Path(base_dir)
    pth_path = Path.home().joinpath('.local', 'lib', 'python3.6', 'site-packages', 'biobb.pth')
    apo = ligand.lower() == 'apo'
    if apo:
        template_yaml_path = base_dir.joinpath('PMX', 'pmxlaunch', 'workflows', 'APO', 'pmx_apo.yaml')
        template_py_path = base_dir.joinpath('PMX', 'pmxlaunch', 'workflows', 'APO', 'pmx_apo.py')
    else:
        template_yaml_path = base_dir.joinpath('PMX', 'pmxlaunch', 'workflows', 'HOLO', 'pmxLig.yaml')
        template_py_path = base_dir.joinpath('PMX', 'pmxlaunch', 'workflows', 'HOLO', 'pmxLig.py')

    # Check if  biobb.pth file exists and if not exists create it
    if not pth_path.exists():
        create_biobb_pth_file(pth_path)

    # Get input trajs
    traj_wt_tpr_path, traj_wt_xtc_path = get_input_tpr_xtc_paths('WT', wt_replica, ligand, base_dir)
    print('WT trajs: ')
    print(traj_wt_tpr_path, traj_wt_xtc_path)
    print('\n')
    traj_mut_tpr_path, traj_mut_xtc_path = get_input_tpr_xtc_paths(mutation, mut_replica, ligand, base_dir)
    print(f'{mutation} trajs: ')
    print(traj_mut_tpr_path, traj_mut_xtc_path)
    print('\n')

    imaged_traj_available = is_imaged_traj_available(traj_wt_xtc_path, traj_mut_xtc_path)

    # Create working dir path
    working_dir_path = base_dir.joinpath('PMX', 'pmxlaunch', 'runs', f"{three_to_one_mutation(mutation)}_{str(mut_replica)}_WT_{str(wt_replica)}", ligand.upper())
    working_dir_path.mkdir(parents=True, exist_ok=True)

    # Check if it's the first launch
    run_number = 0
    run_dir = working_dir_path.joinpath("wf_pmx")
    config_yaml_path = working_dir_path.joinpath(f"pmx.yaml")
    wf_py_path = working_dir_path.joinpath(f"pmx.py")
    launch_path = working_dir_path.joinpath(f"launch.sh")
    job_name = f"pycompss_pmx"
    while run_dir.exists():
        run_number += 1
        run_dir = working_dir_path.joinpath(f"wf_pmx_{str(run_number)}")
        config_yaml_path = working_dir_path.joinpath(f"pmx_{str(run_number)}.yaml")
        wf_py_path = working_dir_path.joinpath(f"pmx_{str(run_number)}.py")
        launch_path = working_dir_path.joinpath(f"launch_{str(run_number)}.sh")
        job_name = f"pycompss_pmx_{str(run_number)}"

    # Copy py file
    shutil.copyfile(template_py_path, wf_py_path)

    # Read yaml template file
    config_dict = get_template_config_dict(template_yaml_path)
    # Update config_dict
    config_dict['working_dir_path'] = str(run_dir)
    mutation_dict = get_mutation_dict(mutation)
    config_dict['mutations']['stateA'] = f"{mutation_dict['wt']}{str(pmx_resnum)}{mutation_dict['mt']}"
    config_dict['mutations']['stateB'] = f"{mutation_dict['mt']}{str(pmx_resnum)}{mutation_dict['wt']}"
    config_dict['input_trajs']['stateA']['input_tpr_path'] = str(traj_wt_tpr_path)
    config_dict['input_trajs']['stateA']['input_traj_path'] = str(traj_wt_xtc_path)
    config_dict['input_trajs']['stateB']['input_tpr_path'] = str(traj_mut_tpr_path)
    config_dict['input_trajs']['stateB']['input_traj_path'] = str(traj_mut_xtc_path)
    config_dict['step1_trjconv']['properties']['skip'] = trjconv_skip
    if apo:
        config_dict['step11_gmx_grompp']['properties']['mdp']['nsteps'] = fe_nsteps
        config_dict['step11_gmx_grompp']['properties']['mdp']['delta-lambda'] =  float(f'{1 / fe_nsteps:.0g}')
    else:
        clean_ligand_name = ligand.upper()
        if '_' in ligand:
            clean_ligand_name = ligand.upper()[:ligand.rfind('_')]
        config_dict['step4_remove_ligand']['properties']['ligand'] = clean_ligand_name
        config_dict['step7_lig_gmx_appendLigand']['paths']['input_itp_path'] = str(base_dir.joinpath('ITPs', f"{clean_ligand_name}.itp"))
        config_dict['step7_lig_gmx_appendLigand']['paths']['input_posres_itp_path'] = str(base_dir.joinpath('ITPs', f"posre_{clean_ligand_name}.itp"))
        config_dict['step7_lig_gmx_appendLigand']['properties']['posres_name'] = f"POSRES_{clean_ligand_name}"
        config_dict['step14_gmx_grompp']['properties']['mdp']['nsteps'] = fe_nsteps

    with open(config_yaml_path, 'w') as config_yaml_file:
        config_yaml_file.write(yaml.dump(config_dict))

    # Create launch
    with open(launch_path, 'w') as launch_file:
        launch_file.write(f"#!/bin/bash\n")
        launch_file.write(f"\n")
        launch_file.write(f"module purge\n")
        launch_file.write(f"\n")
        launch_file.write(f"module load ANACONDA/2018.12_py3\n")
        launch_file.write(f"source activate biobb\n")
        launch_file.write(f"\n")
        launch_file.write(f"# COMPSs environment\n")
        launch_file.write(f"export COMPSS_PYTHON_VERSION=none\n")
        launch_file.write(f"# COMPSs release\n")
        launch_file.write(f"module load COMPSs/{compss_version}\n")
        launch_file.write(f"\n")
        launch_file.write(f"# Singularity\n")
        launch_file.write(f"module load singularity\n")
        launch_file.write(f"\n")
        launch_file.write(f"#GROMACS 2019\n")
        launch_file.write(f"module load intel/2018.4 impi/2018.4 mkl/2018.4 gromacs/2019.1\n")
        launch_file.write(f"\n")
        launch_file.write(f"#Permissions for everyone\n")
        launch_file.write(f"umask ugo+rwx\n")
        launch_file.write(f"\n")
        launch_file.write(f"enqueue_compss ")
        if compss_debug:
            launch_file.write(f"-d ")
        launch_file.write(f"--job_name={job_name} \
                          --num_nodes={num_nodes} \
                          --exec_time={str(time)} \
                          --network=ethernet \
                          --qos={queue}  \
                          {wf_py_path} \
                          --config {config_yaml_path} ")
        if imaged_traj_available:
            launch_file.write(f"--imaged_traj_available ")
        launch_file.write(f"\n")

    subprocess.call(f"bash {launch_path}", shell=True)


def main():
    parser = argparse.ArgumentParser(description="Wrapper of the GROMACS editconf module.")
    parser.add_argument('-m', '--mutation', required=True, help="Mutation in 'Leu858Arg' format")
    parser.add_argument('-prn', '--pmx_resnum', required=False, default=0, type=int, help="(0) [integer]")
    parser.add_argument('-l', '--ligand', required=False, default='APO', type=str, help="(apo) [apo] or [name_of_ligand]")
    parser.add_argument('-wr', '--wt_replica', required=False, default=1, type=int, help="(1) [integer]")
    parser.add_argument('-mr', '--mut_replica', required=False, default=1, type=int, help="(1) [integer]")
    parser.add_argument('-q', '--queue', required=False, default='bsc_ls', type=str, help="(bsc_ls) [bsc_ls|debug]")
    parser.add_argument('-t', '--time', required=False, default=120, type=int, help="(120) [integer] Time in minutes")
    parser.add_argument('-nn', '--num_nodes', required=False, default=1, type=int, help="(1) [integer]")
    parser.add_argument('-cv', '--compss_version', required=False, default='2.6.1', type=str, help="(2.6.1) [version_name]")
    parser.add_argument('-d', '--compss_debug', required=False, help="Compss debug mode", action='store_true')
    parser.add_argument('-fe', '--fe_length', required=False, default=50, type=int, help="(50) [integer] Number of picoseconds")
    parser.add_argument('-nf', '--num_frames', required=False, default=100, type=int, help="(100) [integer] Number of frames to be stracted of trajectory")
    parser.add_argument('--free_energy_dt', required=False, default=0.002, type=float, help="(0.002) [float] Integration time in picoseconds")
    parser.add_argument('--trajectory_total_num_frames', required=False, default=10000, type=int, help="(10000) [integer] Total number of frames of the original trajectory")
    parser.add_argument('--base_dir', required=False, default='/gpfs/projects/bsc23/bsc23513/BioExcel/BioExcel_EGFR_pmx', type=str, help="('/gpfs/projects/bsc23/bsc23513/BioExcel/BioExcel_EGFR_pmx') [path_to_base_dir]")
    args = parser.parse_args()

    # Specific call of each building block
    launch(mutation=args.mutation,
           pmx_resnum=args.pmx_resnum,
           ligand=args.ligand.upper(),
           wt_replica=args.wt_replica,
           mut_replica=args.mut_replica,
           queue=args.queue,
           time=args.time,
           num_nodes=args.num_nodes,
           compss_version=args.compss_version,
           compss_debug=args.compss_debug,
           fe_nsteps=int(args.fe_length/args.free_energy_dt),
           trjconv_skip=args.trajectory_total_num_frames//args.num_frames,
           base_dir=Path(args.base_dir)
           )


if __name__ == '__main__':
    main()
