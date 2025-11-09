import os
import csv
import time
import openvsp as vsp

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)
    return os.path.abspath(path)

def find_parm_id_by_substr(geom_id, substr):
    
    parm_ids = vsp.GetGeomParmIDs(geom_id)
    for pid in parm_ids:
        name = vsp.GetParmName(pid).lower()
        grp = vsp.GetParmGroupName(pid).lower()
        if substr.lower() in name or substr.lower() in grp:
            return pid
    return None

def set_parm_value(geom_id, substr, value):
    
    pid = find_parm_id_by_substr(geom_id, substr)
    if pid:
        vsp.SetParmVal(pid, value)
        return True
    else:
        print(f"WARNING: parm containing '{substr}' not found for geom {geom_id}")
        return False

def list_available_analyses():
    arr = vsp.ListAnalysis()
    print("Available analyses:")
    for a in arr:
        print("  ", a)
    return arr

def find_analysis_by_substr(substr):
    for a in vsp.ListAnalysis():
        if substr.lower() in a.lower():
            return a
    return None

def print_analysis_inputs(analysis):
    print(f"\nInputs for {analysis}:")
    names = vsp.GetAnalysisInputNames(analysis)
    for n in names:
        t = vsp.GetAnalysisInputType(analysis, n)
        print("   ", n, "  (type:", t, ")")
    return names


out_dir = ensure_dir("./out/")
print("Output directory:", out_dir)

vsp.ClearVSPModel()
wing = vsp.AddGeom("WING", "")
print("Created wing:", wing)

root_c = 1
tip_c = 1
span  = 6

set_parm_value(wing, "root", root_c)     
set_parm_value(wing, "tip", tip_c)       
set_parm_value(wing, "span", span)       

set_parm_value(wing, "sweep", 0)
set_parm_value(wing, "dihedral", 0)
set_parm_value(wing, "twist", 0)

vsp.Update()


timestamp = int(time.time())
vsp_filename = os.path.join(out_dir, f"wing_rc{root_c}tc{tip_c}_s{span}{timestamp}.vsp3")
ret = vsp.WriteVSPFile(vsp_filename)
if ret == 0:
    print("Saved VSP file:", vsp_filename)
else:
    print("WriteVSPFile returned non-zero:", ret)