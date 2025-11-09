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
out_dir = os.path.join(out_dir, f"vspaero_{timestamp}")
os.mkdir(out_dir)
vsp_filename = os.path.join(out_dir, f"wing_rc{root_c}tc{tip_c}_s{span}{timestamp}.vsp3")
ret = vsp.WriteVSPFile(vsp_filename)
if ret == 0:
    print("Saved VSP file:", vsp_filename)
else:
    print("WriteVSPFile returned non-zero:", ret)

a_name = find_analysis_by_substr("VSPAEROSweep")
if not a_name:
    print("ERROR: Could not find a VSPAERO analysis by name. Available analyses:")
    list_available_analyses()
    raise SystemExit(1)

print("Using analysis:", a_name)

# Print inputs so you know exact names and types
inputs = print_analysis_inputs(a_name)

# Set defaults then set the inputs we want
vsp.SetAnalysisInputDefaults(a_name)

# We'll try to find and set likely inputs (Mach, Alpha, Reynolds, RefGeom, Iterations)
#  - Mach is usually a double input named 'Mach'
#  - AoA might be AlphaStart/AlphaEnd/AlphaNpts for sweeps OR 'Alpha' for single run
#  - Re can be named 'Re' or 'Reynolds'
#  - Reference geometry input often named 'RefGeom' or 'ReferenceGeom' or 'RefFlag'; it might expect an integer or string
#  - Iterations often named 'NumIter' or 'NIter' or 'Iterations'

def set_analysis_double(analysis, name_substr, value, idx=0):
    names = vsp.GetAnalysisInputNames(analysis)
    for n in names:
        if name_substr.lower() in n.lower():
            v = vsp.GetDoubleAnalysisInput(analysis, n)
            # prepare indata array with same length
            indata = [value] * len(v)
            vsp.SetDoubleAnalysisInput(analysis, n, indata, idx)
            print(f"Set {n} = {indata}")
            return True
    print(f"Warning: analysis double input containing '{name_substr}' not found.")
    return False

def set_analysis_int(analysis, name_substr, intarr, idx=0):
    names = vsp.GetAnalysisInputNames(analysis)
    for n in names:
        if name_substr.lower() in n.lower():
            vsp.SetIntAnalysisInput(analysis, n, intarr, idx)
            print(f"Set {n} = {intarr}")
            return True
    print(f"Warning: analysis int input containing '{name_substr}' not found.")
    return False

def set_analysis_str(analysis, name_substr, strarr, idx=0):
    names = vsp.GetAnalysisInputNames(analysis)
    for n in names:
        if name_substr.lower() in n.lower():
            vsp.SetStringAnalysisInput(analysis, n, strarr, idx)
            print(f"Set {n} = {strarr}")
            return True
    print(f"Warning: analysis string input containing '{name_substr}' not found.")
    return False

# Example: set Mach, single alpha, Re, iterations, reference geometry
set_analysis_double(a_name, "Mach", 0.2)

set_analysis_double(a_name, "AlphaStart", 0)
set_analysis_double(a_name, "AlphaEnd", 10.0)
set_analysis_int(a_name, "AlphaNpts", [10])

set_analysis_double(a_name, "ReCref", 1e7)

# Iterations (search various likely names)
if not set_analysis_int(a_name, "WakeNumIter", [200]):
    if not set_analysis_int(a_name, "WakeNumIter", [200]):
        print("Warning: iterations input not found; solver will use default iterations.")

# Reference geometry: often VSPAERO expects an integer index for a geometry set or string name.
# If the input is a string, set it to the geom name; if it's int, set to the set index.
# We'll try both:
set_analysis_str(a_name, "WingID", [vsp.GetGeomName(wing)])   # tries to set a string ref to our wing name

# OPTIONAL: if you want to run the Compute Geometry step first:
geom_compute_name = find_analysis_by_substr("computegeometry") or find_analysis_by_substr("compute geometry") or "VSPAEROComputeGeometry"
if geom_compute_name:
    # make sure thin geometry set includes our wing; defaults are usually OK
    vsp.SetAnalysisInputDefaults(geom_compute_name)
    vsp.ExecAnalysis(geom_compute_name)
    print("Ran ComputeGeometry step:", geom_compute_name)

# Run the VSPAERO analysis
print("Executing:", a_name)
res_id = vsp.ExecAnalysis(a_name)
print("ExecAnalysis returned results id:", res_id)

# -------------------------
# 3) Extract results and save CSV
# -------------------------
# Find the latest results for VSPAERO (result name varies, inspect results list)
# Use FindLatestResultsID by searching for typical names like 'VSPAERO' in results
all_results = vsp.GetAllResultsNames()
print("All result names available:", all_results)

# Try to find matching results id (there are helper functions: FindLatestResultsID)
# We attempt several common names; if none match we list results and ask to inspect
# candidates = ["VSPAERO", "VSPAEROCase", "VSPAERO_Stab", "VSPAero", "VSPAeroSweep"]
candidates = vsp.GetAllResultsNames()
latest_res_id = None
for cand in candidates:
    try:
        latest_res_id = vsp.FindLatestResultsID(cand)
        if latest_res_id:
            print("Found results for", cand, "->", latest_res_id)
            break
    except Exception:
        pass

if not latest_res_id:
    # Fallback: use the res_id returned by ExecAnalysis
    latest_res_id = res_id
    print("Using ExecAnalysis result id:", latest_res_id)

# Query which data names are present inside this result
if latest_res_id:
    data_names = vsp.GetAllDataNames(latest_res_id)
    print("Data fields in result:", data_names)
    # Common fields: 'CL_tot', 'CD_tot', 'CM', 'CL', 'CD' etc (names vary by version)
    desired = ["CL", "CD", "CL_tot", "CD_tot", "Cm", "Cmz", "L", "D"]
    row = {"case": os.path.basename(vsp_filename)}
    for field in data_names:
        # get double results where applicable
        try:
            dtype = vsp.GetResultsType(latest_res_id, field)
        except Exception:
            dtype = None
        # If it's double(s), extract
        try:
            vals = vsp.GetDoubleResults(latest_res_id, field)
            # store as comma-joined string or single value
            row[field] = ",".join(map(str, vals)) if len(vals) > 1 else (vals[0] if vals else "")
        except Exception:
            # try string results
            try:
                sval = vsp.GetStringResults(latest_res_id, field)
                row[field] = ",".join(sval)
            except Exception:
                row[field] = ""
    # Write a CSV summary (one-row)
    csv_file = os.path.join(out_dir, f"vspaero_results_{timestamp}.csv")
    with open(csv_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["case"] + list(row.keys()))
        # write header
        writer.writeheader()
        writer.writerow(row)
    print("Saved results CSV:", csv_file)
else:
    print("No results id found. Use vsp.PrintResultsDoc or vsp.PrintResults to debug.")

print("ALL DONE.")