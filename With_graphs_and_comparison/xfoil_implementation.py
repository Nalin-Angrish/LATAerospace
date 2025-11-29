import subprocess
import numpy as np
import os

xfoil_path = r"C:\Users\ISHAN\Downloads\XFOIL6.99\xfoil.exe"
def run_xfoil(airfoil, alpha_range, Re=1e6, Mach=0.0, n_iter=100, xfoil_path=xfoil_path):
    # Make a temp command file for XFOIL
    cmd_file = "xfoil_cmd.in"
    polar_file = "polar.txt"
    
    with open(cmd_file, 'w') as f:
        # Load airfoil
        if airfoil.lower().startswith("naca"):
            f.write(f"NACA {airfoil[4:]}\n")
        else:
            f.write(f"LOAD {airfoil}\n")
            f.write(f"{airfoil}\n")  # confirm name
        
        # Oper menu
        f.write("OPER\n")
        f.write(f"VISC {Re}\n")
        f.write(f"MACH {Mach}\n")
        f.write(f"ITER {n_iter}\n")
        
        # Set polar output
        f.write(f"PACC\n{polar_file}\n\n")
        
        # Sweep AoA
        for alpha in alpha_range:
            f.write(f"ALFA {alpha}\n")
        
        # Close polar file
        f.write("PACC\n")
        f.write("QUIT\n")
    
    # Run XFOIL
    p = subprocess.Popen(xfoil_path,
                         stdin=open(cmd_file, 'r'),
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         shell=True)
    out, err = p.communicate()

    # Read results
    data = np.loadtxt(polar_file, skiprows=12)  # first 12 lines are header
    return data


# Example usage:
if __name__ == "__main__":
    alphas = np.linspace(-4, 12, 17)
    polar = run_xfoil("NACA2412", alphas)
    print("alpha, Cl, Cd, Cm:")
    print(polar[:, :4])  # first 4 columns
