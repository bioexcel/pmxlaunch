import props

def main(cfgs):
    from wfs import pmx_wf
    for cfg in cfgs:
       pmx_wf("--config", cfg, "--imaged_traj_available")

if __name__ == '__main__':
    #read properties and set
    props.set_properties("2","48","/gpfs/scratch/bsc19/bsc19611/test_bioExcel/","-d")
    cfgs = ["/gpfs/scratch/bsc19/bsc19611/test_bioExcel/IRE/pmx.yaml", "/gpfs/scratch/bsc19/bsc19611/test_bioExcel/IRE2/pmx.yaml"]   
    main(cfgs)

