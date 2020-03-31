import props

def main(cfgs):
    from wfs import pmx_wf
    for cfg in cfgs:
       pmx_wf("--config", cfg, "--imaged_traj_available")

if __name__ == '__main__':
    #read properties and set
    props.set_properties("2","48","/gpfs/projects/bsc23/bsc23513/BioExcel/BioExcel_EGFR_pmx/PMX/pmxlaunch/runs/G796S_1_WT_1_10f_50ps/alltoall_alltoall/IRE/","-d")
    cfgs = ["/gpfs/projects/bsc23/bsc23513/BioExcel/BioExcel_EGFR_pmx/PMX/pmxlaunch/runs/G796S_1_WT_1_10f_50ps/alltoall_alltoall/IRE/pmx.yaml"]   
    main(cfgs)

