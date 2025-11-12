"""
Run a VSPAero analysis on a simple wing geometry created in OpenVSP.
Use this script to optimize wing parameters c and b over specified ranges.
"""
import os
import numpy as np
import pandas as pd
import openvsp as vsp

## Define the wing geometry parameter ranges:
range_root_c = np.linspace(0.5, 2.0, 4)  # Root chord from 0.5 to 2.0
range_span   = np.linspace(1.0, 5.0, 20)
sample_space = [(c, b) for c in range_root_c for b in range_span]


# For now, we're using the NACA 2412 airfoil.
# In the future, we may integrate airfoil optimization directly.
TC=0.12
Camber=0.01562
Camber_Location=0.40

def ensure_dir(path):
    """
    Ensure that a directory exists. Create it if it doesn't.
    """
    os.makedirs(path, exist_ok=True)
    return os.path.abspath(path)

out_dir = ensure_dir("../out/openvsp_simplewing_optimizer")
csv_dir = ensure_dir(os.path.join(out_dir, "results"))
print("Output directory:", out_dir)

for point in sample_space:
    ## Define the wing geometry parameters
    root_c = point[0]
    tip_c = root_c # Assuming a rectangular wing
    span  = point[1]

    vsp.VSPRenew()
    vsp.ClearVSPModel()

    wing = vsp.AddGeom("WING", "")
    vsp.SetGeomName(wing, "MyWing")

    print("Created wing:", wing)
    vsp.Update()

    # Set the wing geometry parameters
    vsp.SetParmVal(wing, "Root_Chord", "XSec_1", root_c)
    vsp.SetParmVal(wing, "Tip_Chord", "XSec_1", tip_c)
    vsp.SetParmVal(wing, "Span", "XSec_1", span)
    vsp.SetParmVal(wing, "Sweep", "XSec_1", 0)

    vsp.Update()

    xsec_surf = vsp.GetXSecSurf(wing, 0)

    vsp.ChangeXSecShape(xsec_surf, 0, vsp.XS_FOUR_SERIES)
    vsp.Update()
    root_xsec = vsp.GetXSec(xsec_surf, 0)
    thick_parm = vsp.GetXSecParm(root_xsec, "ThickChord")
    camber_parm = vsp.GetXSecParm(root_xsec, "Camber")
    camber_loc_parm = vsp.GetXSecParm(root_xsec, "CamberLoc")

    vsp.SetParmVal(thick_parm, TC)
    vsp.SetParmVal(camber_parm, Camber)
    vsp.SetParmVal(camber_loc_parm, Camber_Location)
    vsp.ChangeXSecShape(xsec_surf, 1, vsp.XS_FOUR_SERIES)
    vsp.Update()

    tip_xsec = vsp.GetXSec(xsec_surf, 1)
    thick_parm = vsp.GetXSecParm(tip_xsec, "ThickChord")
    camber_parm = vsp.GetXSecParm(tip_xsec, "Camber")
    camber_loc_parm = vsp.GetXSecParm(tip_xsec, "CamberLoc")
    vsp.SetParmVal(thick_parm, TC)
    vsp.SetParmVal(camber_parm, Camber)
    vsp.SetParmVal(camber_loc_parm, Camber_Location)
    vsp.Update()

    vsp.WriteVSPFile("../out/openvsp_simplewing_optimizer/simple_wing.vsp3")

    # VSPAero Setup

    vsp.SetParmVal(wing, "Sym_Planar_Flag", "Sym", 1)
    vsp.SetParmVal(wing, "Sym_Ancestor", "Sym", 0)

    analysis_name = "VSPAEROComputeGeometry"
    vsp.SetAnalysisInputDefaults(analysis_name)
    vsp.SetIntAnalysisInput(analysis_name, "GeomSet", [0], 0)

    print("\nGenerating VSPAero geometry...")
    comp_geom_results = vsp.ExecAnalysis(analysis_name)
    print("Geometry generation complete")

    vspaero_analysis = "VSPAEROSweep"
    vsp.SetAnalysisInputDefaults(vspaero_analysis)


    ref_area = 2*root_c * span
    ref_chord = root_c
    ref_span = 2*span
    mach=0.06
    n=15 # Here change the number of iterations to converge the wake

    vsp.SetDoubleAnalysisInput(vspaero_analysis, "Sref", [ref_area], 0)
    vsp.SetDoubleAnalysisInput(vspaero_analysis, "bref", [ref_span], 0)
    vsp.SetDoubleAnalysisInput(vspaero_analysis, "cref", [ref_chord], 0)

    vsp.SetDoubleAnalysisInput(vspaero_analysis, "Xcg", [0], 0)
    vsp.SetDoubleAnalysisInput(vspaero_analysis, "Ycg", [0.0], 0)
    vsp.SetDoubleAnalysisInput(vspaero_analysis, "Zcg", [0.0], 0)

    # Here we set the angle of attack and sideslip ranges and number of points
    vsp.SetDoubleAnalysisInput(vspaero_analysis, "AlphaStart", [0.0], 0)
    vsp.SetDoubleAnalysisInput(vspaero_analysis, "AlphaEnd", [10.0], 0)
    vsp.SetIntAnalysisInput(vspaero_analysis, "AlphaNpts", [11], 0)

    vsp.SetDoubleAnalysisInput(vspaero_analysis, "BetaStart", [0.0], 0)
    vsp.SetDoubleAnalysisInput(vspaero_analysis, "BetaEnd", [0.0], 0)
    vsp.SetIntAnalysisInput(vspaero_analysis, "BetaNpts", [1], 0)

    vsp.SetDoubleAnalysisInput(vspaero_analysis, "MachStart", [mach], 0)
    vsp.SetDoubleAnalysisInput(vspaero_analysis, "MachEnd", [0], 0)
    vsp.SetIntAnalysisInput(vspaero_analysis, "MachNpts", [1], 0)
    print("Using alternative parameter names...")

    vsp.SetIntAnalysisInput(vspaero_analysis, "WakeNumIter", [n], 0)

    results_id = vsp.ExecAnalysis(vspaero_analysis)
    print("VSPAero analysis complete!")
    #--------------------------------------------------------------------------

    res_id = vsp.FindLatestResultsID("VSPAERO_Polar")
    if not res_id:
        res_id = results_id
    if res_id:
        data_names = vsp.GetAllDataNames(res_id)
        print("Data fields in result:", data_names)
        alpha = vsp.GetDoubleResults(res_id, "Alpha", 0)
        CLtot = vsp.GetDoubleResults(res_id, "CLtot", 0)
        CDtot = vsp.GetDoubleResults(res_id, "CDtot", 0)
        L_D = vsp.GetDoubleResults(res_id, "L_D", 0)
        path = os.path.join(csv_dir, f"simple_wing_{root_c}_{span}.csv")
        pd.DataFrame({
            "Alpha": alpha,
            "CLtot": CLtot,
            "CDtot": CDtot,
            "L/D": L_D
        }).to_csv(path, index=False)
        print(f"Results saved to {path}")
    else:
        print("No results id found. Use vsp.PrintResultsDoc or vsp.PrintResults to debug.")
