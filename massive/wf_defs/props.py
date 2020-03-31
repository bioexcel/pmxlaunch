cn="1"
cu="24"
wf_path="./"
compss_flags="-d"

def set_properties(comp_nodes, comp_units, path, flags):
    set_cn(comp_nodes)
    set_cu(comp_units)
    set_wf_path(path)
    set_flags(flags)

def set_cn(comp_nodes):
    global cn
    cn=comp_nodes

def get_cn():
    global cn
    return cn

def get_cu():
    global cu
    return cu

def set_cu(comp_units):
    global cu
    cu=comp_units

def set_wf_path(path):
    global wf_path
    wf_path = path

def get_wf_path():
    global wf_path
    return wf_path

def set_flags(flags):
    global compss_flags
    compss_flags = flags

def get_flags():
    global compss_flags
    return compss_flags

