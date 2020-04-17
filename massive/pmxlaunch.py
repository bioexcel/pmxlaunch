#!/usr/bin/env python3

from pathlib import Path
import argparse
from Bio.PDB.Polypeptide import three_to_one
from Bio.PDB.Polypeptide import one_to_three
import re
import shutil
import oyaml as yaml
from csv_util import parse_csv
from wf_defs import props


def create_biobb_pth_file(file_path):
    with open(file_path, 'w') as pth_file:
        pth_file.write("/home/bsc23/bsc23210/macshare/biobb_common \n")
        pth_file.write("/home/bsc23/bsc23210/macshare/biobb_md \n")
        pth_file.write("/home/bsc23/bsc23210/macshare/biobb_pmx \n")
        pth_file.write("/home/bsc23/bsc23210/macshare/biobb_analysis \n")
        pth_file.write("/home/bsc23/bsc23210/macshare/biobb_adapters \n")
        pth_file.write("/home/bsc23/bsc23210/macshare/biobb_wf_pmxligand \n")
        pth_file.write("/home/bsc23/bsc23210/macshare/biobb_structure_utils")


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

def apply_defaults(conf_dict):
    defaults_dict = {
        'pmx_resnum': int(get_mutation_dict(conf_dict['mutation'])['resnum']) - 693,
        'ligand': 'APO',
        'wt_replica': 1,
        'mut_replica': 1,
        'fe_length': 50,
        'num_frames': 100,
        'free_energy_dt': 0.002,
        'trajectory_total_num_frames': 10000,
        'base_dir': '/gpfs/projects/bsc23/bsc23513/BioExcel/BioExcel_EGFR_pmx'
    }
    for key, value in conf_dict:
        if not value:
            conf_dict[key] = defaults_dict[key]


def launch(mutation, pmx_resnum, ligand, wt_replica, mut_replica, fe_nsteps, fe_length,
           base_dir, output_dir, fe_dt, num_frames, wt_trjconv_skip, mut_trjconv_skip,
           wt_start, wt_end, mut_start, mut_end):


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
    long_name = f"{three_to_one_mutation(mutation)}_{str(mut_replica)}_WT_{str(wt_replica)}_{str(num_frames)}f_{str(fe_length)}ps"
    wt_start_str = str(wt_start)
    wt_end_str = str(wt_end)+'ps'
    if wt_start_str == '0' and wt_end_str == '0ps':
        wt_start_str = 'all'
        wt_end_str = 'all'
    mut_start_str = str(mut_start)
    mut_end_str = str(mut_end) + 'ps'
    if mut_start_str == '0' and mut_end_str == '0ps':
        mut_start_str = 'all'
        mut_end_str = 'all'
    working_dir_path = base_dir.joinpath('PMX', 'pmxlaunch', 'runs', long_name, f"{wt_start_str}to{wt_end_str}_{mut_start_str}to{mut_end_str}", ligand.upper())

    if output_dir:
        if output_dir.startswith('/'):
            working_dir_path = Path(output_dir).resolve()
        else:
            working_dir_path = base_dir.joinpath('PMX', 'pmxlaunch', 'runs', output_dir)
    working_dir_path.mkdir(parents=True, exist_ok=True)

    # Check if it's the first launch
    run_number = 0
    run_dir = working_dir_path.joinpath("wf_pmx")
    config_yaml_path = working_dir_path.joinpath(f"pmx.yaml")
    wf_py_path = working_dir_path.joinpath(f"pmx.py")
    while run_dir.exists():
        run_number += 1
        run_dir = working_dir_path.joinpath(f"wf_pmx_{str(run_number)}")
        config_yaml_path = working_dir_path.joinpath(f"pmx_{str(run_number)}.yaml")
        wf_py_path = working_dir_path.joinpath(f"pmx_{str(run_number)}.py")

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
    config_dict['step1_trjconv_stateA']['properties']['skip'] = wt_trjconv_skip
    config_dict['step1_trjconv_stateA']['properties']['start'] = wt_start
    config_dict['step1_trjconv_stateA']['properties']['end'] = wt_end
    config_dict['step1_trjconv_stateB']['properties']['skip'] = mut_trjconv_skip
    config_dict['step1_trjconv_stateB']['properties']['start'] = mut_start
    config_dict['step1_trjconv_stateB']['properties']['end'] = mut_end

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
        config_dict['step14_gmx_grompp']['properties']['mdp']['dt'] = fe_dt

    with open(config_yaml_path, 'w') as config_yaml_file:
        config_yaml_file.write(yaml.dump(config_dict))

    from wf_defs.wfs import lig_wf, apo_wf
    wf = lig_wf
    if apo:
        wf = apo_wf

    if imaged_traj_available:
        wf("--config", config_yaml_path, "--imaged_traj_available")
    else:
        wf("--config", config_yaml_path,)


def main():
    parser = argparse.ArgumentParser(description="Create configuration files and choose which workflow should be launched")
    parser.add_argument('--csv', required=True, help="Config CSV file")
    # parser.add_argument('-m', '--mutation', required=True, help="Mutation in 'Leu858Arg' format")
    # parser.add_argument('-prn', '--pmx_resnum', required=False, default=0, type=int, help="(0) [integer]")
    # parser.add_argument('-l', '--ligand', required=False, default='APO', type=str, help="(apo) [apo] or [name_of_ligand]")
    # parser.add_argument('-wr', '--wt_replica', required=False, default=1, type=int, help="(1) [integer]")
    # parser.add_argument('-mr', '--mut_replica', required=False, default=1, type=int, help="(1) [integer]")
    # parser.add_argument('-q', '--queue', required=False, default='bsc_ls', type=str, help="(bsc_ls) [bsc_ls|debug]")
    # parser.add_argument('-t', '--time', required=False, default=120, type=int, help="(120) [integer] Time in minutes")
    # parser.add_argument('-nn', '--num_nodes', required=False, default=1, type=int, help="(1) [integer]")
    # parser.add_argument('-cv', '--compss_version', required=False, default='2.6.1', type=str, help="(2.6.1) [version_name]")
    # parser.add_argument('-d', '--compss_debug', required=False, help="Compss debug mode", action='store_true')
    # parser.add_argument('-fe', '--fe_length', required=False, default=50, type=int, help="(50) [integer] Number of picoseconds")
    # parser.add_argument('-nf', '--num_frames', required=False, default=100, type=int, help="(100) [integer] Number of frames to be stracted of trajectory")
    # parser.add_argument('--free_energy_dt', required=False, default=0.002, type=float, help="(0.002) [float] Integration time in picoseconds")
    # parser.add_argument('--trajectory_total_num_frames', required=False, default=10000, type=int, help="(10000) [integer] Total number of frames of the original trajectory")
    # parser.add_argument('--base_dir', required=False, default='/gpfs/projects/bsc23/bsc23513/BioExcel/BioExcel_EGFR_pmx', type=str, help="('/gpfs/projects/bsc23/bsc23513/BioExcel/BioExcel_EGFR_pmx') [path_to_base_dir]")
    args = parser.parse_args()

    # Read CSV to list of dicts
    execution_dict_list, field_names = parse_csv(args.csv)
    print('Dict List of CSV:')
    print(execution_dict_list)
    print('\n\n\n\n\n\n\n\n')

    print('Field names of CSV:')
    print(field_names)
    print('\n\n\n\n\n\n\n\n')

    # Properties of each subworkflow_launch
    props.set_properties("2", "48", "/gpfs/projects/bsc23/bsc23513/BioExcel/BioExcel_EGFR_pmx/PMX/pmxlaunch/workflows/", "-d")

    for line_num, execution_conf in enumerate(execution_dict_list):
        # Apply defaults and launch
        print(f'Apply defaults to configuration dict {line_num}:')
        print(execution_conf)
        print('\n\n\n\n\n\n\n\n')

        print(f'Launch function {line_num}:')
        print(f"launch(mutation={execution_conf['mutation']},\
               pmx_resnum={execution_conf['pmx_resnum']}, \
               ligand={execution_conf['ligand']}, \
               wt_replica={execution_conf['wt_replica']}, \
               mut_replica={execution_conf['mut_replica']}, \
               fe_nsteps={int(execution_conf['fe_length']/execution_conf['free_energy_dt'])}, \
               trjconv_skip={execution_conf['trajectory_total_num_frames']//execution_conf['num_frames']}, \
               base_dir={Path(execution_conf['base_dir'])} \
               )")
        print('\n\n\n\n\n\n\n\n')

        launch(mutation=execution_conf['mutation'],
               pmx_resnum=execution_conf['pmx_resnum'],
               ligand=execution_conf['ligand'],
               wt_replica=execution_conf['wt_replica'],
               mut_replica=execution_conf['mut_replica'],
               fe_nsteps=int(execution_conf['fe_length']/execution_conf['free_energy_dt']),
               fe_length=execution_conf['fe_length'],
               output_dir=execution_conf['output_dir'],
               fe_dt=execution_conf['free_energy_dt'],
               num_frames=execution_conf['num_frames'],
               wt_trjconv_skip=execution_conf['wt_start_end_num_frames'] // execution_conf['num_frames'],
               mut_trjconv_skip=execution_conf['mut_start_end_num_frames'] // execution_conf['num_frames'],
               wt_start=execution_conf['wt_start'],
               wt_end=execution_conf['wt_end'],
               mut_start=execution_conf['mut_start'],
               mut_end=execution_conf['mut_end'],
               base_dir=Path(execution_conf['base_dir'])
               )


if __name__ == '__main__':
    main()
