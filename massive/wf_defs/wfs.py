from pycompss.api.task import task
from pycompss.api.compss import compss
from pycompss.api.constraint import constraint
from pycompss.api.parameter import *
import props

@compss(flags=props.get_flags(), app_name=props.get_wf_path()+"HOLO/pmxLig.py", working_dir="$PWD", computing_nodes=props.get_cn(), fail_by_exit_value=True)
@constraint(computing_units=props.get_cu())
@task(config_file=FILE_IN, on_failure="IGNORE")
def lig_wf(config_flag, config_file, imaged_trajectories_flag):
    pass

@compss(flags=props.get_flags(), app_name=props.get_wf_path()+"APO/pmx_apo.py", computing_nodes=props.get_cn())
@constraint(computing_units=props.get_cu())
@task(config_file=FILE_IN, on_failure="IGNORE")
def apo_wf(config_flag, config_file, imaged_trajectories_flag):
    pass

