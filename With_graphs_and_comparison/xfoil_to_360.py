#!/usr/bin/env python3
"""
xfoil_to_360.py

Usage:
    - Place xfoil executable somewhere and set xfoil_path accordingly (default 'xfoil' / 'xfoil.exe').
    - Call main() or import XfoilPolar360 and use programmatically.

What it does:
    1. Launches XFOIL with a command script (viscous polar sweep).
    2. Parses the produced polar file.
    3. Extrapolates to 0-360° using Viterna-Corrigan.
    4. Plots Cl, Cd, Cm vs alpha and saves a CSV for BEMT use.
"""

import subprocess
import numpy as np
import pandas as pd
import os
import tempfile
import sys
import matplotlib.pyplot as plt
from typing import Tuple, Optional, Sequence

# -------------------------
# Utility: run XFOIL
# -------------------------
xfoil_path = r"C:\Users\ISHAN\Downloads\XFOIL6.99\xfoil.exe"
def run_xfoil_and_get_polar(airfoil: str,
                            alpha_range: Sequence[float],
                            Re: float = 1e6,
                            Mach: float = 0.0,
                            n_iter: int = 200,
                            xfoil_path: str = xfoil_path,
                            polar_filename: str = "polar.txt",
                            verbose: bool = True) -> str:
    """
    Run XFOIL with a generated command file and return the path to polar file.
    airfoil: "NACA2412" or path to .dat file (e.g. "mh114.dat")
    alpha_range: iterable of AoA in degrees (e.g. np.linspace(-5, 15, 21))
    Re: Reynolds number
    xfoil_path: executable path (e.g. "C:\\xfoil\\xfoil.exe" on Windows)
    polar_filename: name of output polar file that XFOIL will write
    Returns: path to polar file (polar_filename)
    """
    # Create a temporary command file for XFOIL
    tmp_cmd = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.in')
    cmd_path = tmp_cmd.name

    # Build commands
    def write(s):
        tmp_cmd.write(s + "\n")

    # Load airfoil
    if airfoil.strip().lower().startswith("naca"):
        # example: "NACA2412"
        code = airfoil.strip()[4:].strip()
        write(f"NACA {code}")
    else:
        # assume a file path
        # supply full path in LOAD then the name for confirmation
        fullp = os.path.abspath(airfoil)
        write(f"LOAD {fullp}")
        # XFOIL expects the airfoil identifier on the next line (some versions)
        write(os.path.basename(fullp))

    write("PLOP")           # disable plotting to console if possible
    write("G")              # send 'g' to disable graphics (some versions)
    write("")
    write("OPER")
    # viscous mode with Reynolds number
    write(f"VISC {Re:g}")
    write(f"MACH {Mach:g}")
    write(f"ITER {n_iter}")
    # Turn off dump of every iteration
    write("PACC")
    write(polar_filename)
    write("")  # second line blank as required by XFOIL PACC

    # Sweep alphas
    for a in alpha_range:
        write(f"ALFA {float(a):.6g}")

    # Close PACC and quit
    write("PACC")
    write("QUIT")

    tmp_cmd.flush()
    tmp_cmd.close()

    if verbose:
        print(f"[xfoil] Running {xfoil_path} with command file {cmd_path} ...")
    # Launch XFOIL
    # On Windows, setting shell=True often helps if path includes backslashes and exe is direct.
    shell_flag = (os.name == 'nt')

    try:
        with open(cmd_path, 'r') as cmdf:
            # We give the command file as stdin for xfoil
            p = subprocess.Popen([xfoil_path], stdin=cmdf, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=shell_flag)
            out, err = p.communicate(timeout=60 + len(alpha_range)*2)
    except FileNotFoundError:
        raise RuntimeError(f"XFOIL executable not found at '{xfoil_path}'. Provide correct path.")
    except subprocess.TimeoutExpired:
        p.kill()
        out, err = p.communicate()
        raise RuntimeError("XFOIL run timed out.")

    if verbose:
        stdout = out.decode(errors='ignore') if isinstance(out, (bytes, bytearray)) else str(out)
        stderr = err.decode(errors='ignore') if isinstance(err, (bytes, bytearray)) else str(err)
        print("[xfoil] stdout snippet:")
        print("\n".join(stdout.splitlines()[:15]))
        if stderr.strip():
            print("[xfoil] stderr:")
            print(stderr)

    # cleanup command file
    try:
        os.remove(cmd_path)
    except Exception:
        pass

    if not os.path.exists(polar_filename):
        raise RuntimeError(f"Expected polar file '{polar_filename}' not created by XFOIL.")

    return polar_filename

# -------------------------
# Utility: parse XFOIL polar robustly
# -------------------------
def parse_xfoil_polar(polar_path: str) -> pd.DataFrame:
    """
    Parse an XFOIL polar file into a pandas DataFrame.
    Attempts to find the header (line with 'alpha' or '-----') and read from there.
    Typical columns: alpha CL CD CDp CM Top_Xtr Bot_Xtr
    """
    with open(polar_path, 'r', errors='ignore') as f:
        lines = f.readlines()

    # find header line index (a line that contains 'alpha' or '-----' separator)
    header_idx = None
    for i, line in enumerate(lines):
        low = line.lower()
        if 'alpha' in low and 'cl' in low:
            header_idx = i
            break

    if header_idx is None:
        # fallback: find the first line that looks numeric (two+ columns)
        for i, line in enumerate(lines):
            parts = line.strip().split()
            if len(parts) >= 3:
                try:
                    [float(x) for x in parts[:3]]
                    header_idx = i - 1  # assume header is previous line
                    break
                except Exception:
                    continue

    if header_idx is None:
        # as last resort, attempt to read skipping 12 rows (common)
        skip = 12
        try:
            df = pd.read_csv(polar_path, delim_whitespace=True, header=None, skiprows=skip)
            df.columns = ['alpha', 'CL', 'CD', 'CDp', 'CM', 'Top_Xtr', 'Bot_Xtr'][:df.shape[1]]
            return df
        except Exception as e:
            raise RuntimeError("Could not find header in polar file and fallback failed.") from e

    # header line content -> derive column names
    header_line = lines[header_idx].strip()
    colnames = header_line.split()
    # rows after header until blank or until end
    data_lines = []
    for line in lines[header_idx+1:]:
        if not line.strip():
            break
        # If line begins with non-digit and non-minus maybe it's a footer
        parts = line.strip().split()
        # lines with numeric values -> keep
        try:
            float(parts[0])
            data_lines.append(line)
        except Exception:
            # stop if not numeric
            break

    # create a temporary string to read via pandas
    from io import StringIO
    data_str = "".join(data_lines)
    if not data_str.strip():
        raise RuntimeError("No numeric data found in polar file after header.")
    df = pd.read_csv(StringIO(data_str), delim_whitespace=True, header=None)
    # assign column names (if sizes mismatch, create generic names)
    if len(colnames) >= df.shape[1]:
        df.columns = colnames[:df.shape[1]]
    else:
        # create generic column names
        df.columns = ['alpha', 'CL', 'CD', 'CDp', 'CM', 'Top_Xtr', 'Bot_Xtr'][:df.shape[1]]
    return df

# -------------------------
# Viterna-Corrigan extrapolation
# -------------------------
def viterna_corrigan_extrapolate(alpha_deg: np.ndarray,
                                 CL_xfoil: np.ndarray,
                                 CD_xfoil: np.ndarray,
                                 AR: float = 6.0,
                                 alpha_stall_pos: Optional[float] = None,
                                 alpha_stall_neg: Optional[float] = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Extrapolate to 0-360 deg using Viterna-Corrigan.
    alpha_deg: array of angles (deg) where XFOIL data exists (typically -10..+15)
    CL_xfoil, CD_xfoil: corresponding arrays
    AR: aspect ratio for CD_max estimate (use a large AR for thin infinite wing behaviour)
    alpha_stall_pos / alpha_stall_neg: optionally provide stall angles; else estimate from max(CL).
    Returns:
        alpha_full (0..360), CL_full, CD_full
    """
    # 1) Determine CL_max (positive stall) and stall alpha
    idx_max = np.nanargmax(CL_xfoil)
    CL_max = CL_xfoil[idx_max]
    alpha_CLmax = alpha_deg[idx_max]

    # estimate stall angles if not provided
    if alpha_stall_pos is None:
        alpha_stall_pos = alpha_CLmax
    if alpha_stall_neg is None:
        idx_min = np.nanargmin(CL_xfoil)
        alpha_stall_neg = alpha_deg[idx_min]

    # Viterna parameters
    # CD_max for a rough flat-plate style
    CD_max = (1.11 + 0.018 * AR) / 2.0

    # Prepare full alpha array 0..360
    alpha_full = np.arange(0, 361)
    CL_full = np.zeros_like(alpha_full, dtype=float)
    CD_full = np.zeros_like(alpha_full, dtype=float)

    # Put the measured XFOIL region into the full arrays (for the exact alpha points)
    # Round XFOIL alphas to nearest integers to place them; but also keep an interpolation later
    # We'll interpolate XFOIL CL/CD across -180..180 domain then overwrite central region
    # Map XFOIL angles into -180..180
    a_x = np.array(alpha_deg, dtype=float)
    # normalize to -180..180
    a_x_norm = ((a_x + 180) % 360) - 180

    # create fine interpolation over -180..180 from XFOIL data
    # For interpolation we need sorted arrays
    sort_idx = np.argsort(a_x_norm)
    a_sorted = a_x_norm[sort_idx]
    CL_sorted = CL_xfoil[sort_idx]
    CD_sorted = CD_xfoil[sort_idx]

    # clip to avoid duplicates
    # build interpolation functions with linear extrapolation allowed
    from scipy.interpolate import interp1d
    interp_CL = interp1d(a_sorted, CL_sorted, kind='linear', fill_value='extrapolate', assume_sorted=True)
    interp_CD = interp1d(a_sorted, CD_sorted, kind='linear', fill_value='extrapolate', assume_sorted=True)

    # We will treat the "attached" region as near the XFOIL sweep: find min/max of a_sorted
    a_min, a_max = a_sorted.min(), a_sorted.max()

    for i, a in enumerate(alpha_full):
        # map a to -180..180
        a_norm = ((a + 180) % 360) - 180

        # If a_norm is within the XFOIL region +/- a margin (like +/-2 deg), use interpolated XFOIL
        margin = 2.0
        if (a_norm >= a_min - margin) and (a_norm <= a_max + margin):
            CL_full[i] = float(interp_CL(a_norm))
            CD_full[i] = float(interp_CD(a_norm))
        else:
            # Viterna-Corrigan formulas (post-stall / separated flow)
            alpha_rad = np.deg2rad(a)
            # CL: CL_max * 0.5 * sin(2 alpha)
            CL_vit = CL_max * 0.5 * np.sin(2.0 * alpha_rad)
            # CD: CD_max * sin(alpha)^2
            CD_vit = CD_max * (np.sin(alpha_rad) ** 2)
            CL_full[i] = CL_vit
            CD_full[i] = CD_vit

    # Ensure continuity near 0/360 by smoothing small overlaps (simple moving average on boundaries)
    # optional smoothing:
    def smooth_edge(arr):
        arr[:3] = np.convolve(arr[:6], np.ones(6)/6.0, mode='same')[:3]
        arr[-3:] = np.convolve(arr[-6:], np.ones(6)/6.0, mode='same')[-3:]
        return arr

    CL_full = smooth_edge(CL_full)
    CD_full = smooth_edge(CD_full)

    return alpha_full, CL_full, CD_full

# -------------------------
# Helper: generate Cm by simple extrapolation if not present
# -------------------------
def estimate_Cm_full(alpha_full: np.ndarray, df_polar: pd.DataFrame) -> np.ndarray:
    """
    If Cm exists in XFOIL polar, we interpolate it; otherwise produce a simple synthetic Cm curve:
    small linear slope near zero and decays to near zero in fully separated region.
    """
    if 'CM' in df_polar.columns:
        # interpolate CM over -180..180 domain
        a_src = df_polar['alpha'].values
        cm_src = df_polar['CM'].values
        # map to -180..180
        a_src_norm = ((a_src + 180) % 360) - 180
        from scipy.interpolate import interp1d
        fcm = interp1d(a_src_norm, cm_src, kind='linear', fill_value='extrapolate', assume_sorted=False)
        cm_full = np.array([float(fcm(((a + 180) % 360) - 180)) for a in alpha_full])
        return cm_full
    else:
        # synthetic Cm: start at small negative value and reduce magnitude away from small alphas
        cm0 = -0.02  # typical small camber moment
        cm_full = cm0 * np.exp(- (np.sin(np.deg2rad(alpha_full))**2) * 4.0)
        return cm_full

# -------------------------
# Main class wrapper
# -------------------------
class XfoilPolar360:
    def __init__(self, xfoil_path: str = "xfoil", verbose: bool = True):
        self.xfoil_path = xfoil_path
        self.verbose = verbose

    def make_polar_and_extrapolate(self,
                                   airfoil: str,
                                   alphas: Sequence[float],
                                   Re: float = 1e6,
                                   Mach: float = 0.0,
                                   AR: float = 6.0,
                                   polar_filename: str = "polar.txt",
                                   out_csv: str = "polar_0_360.csv",
                                   plot_fig: str = "polar_0_360.png"):
        # 1) Run XFOIL
        polar_path = run_xfoil_and_get_polar(airfoil=airfoil,
                                             alpha_range=alphas,
                                             Re=Re,
                                             Mach=Mach,
                                             xfoil_path=self.xfoil_path,
                                             polar_filename=polar_filename,
                                             verbose=self.verbose)
        # 2) Parse polar
        df = parse_xfoil_polar(polar_path)
        if self.verbose:
            print(f"[parse] Found columns: {list(df.columns)}")
            print(df.head())

        # ensure alpha column exists
        if 'alpha' not in [c.lower() for c in df.columns]:
            # try standard names
            df.columns = [c.strip() for c in df.columns]
        # normalize column names to lowercase for access
        cols = {c.lower(): c for c in df.columns}
        alpha_col = cols.get('alpha', None)
        cl_col = cols.get('cl', None)
        cd_col = cols.get('cd', None)
        if alpha_col is None or cl_col is None or cd_col is None:
            raise RuntimeError("Parsed polar does not contain alpha/CL/CD columns.")

        alpha = df[alpha_col].values.astype(float)
        CL = df[cl_col].values.astype(float)
        CD = df[cd_col].values.astype(float)

        # 3) Extrapolate
        # We rely on scipy for interpolation inside the viterna function
        try:
            alpha_full, CL_full, CD_full = viterna_corrigan_extrapolate(alpha, CL, CD, AR=AR)
        except Exception as e:
            raise RuntimeError("Extrapolation failed; ensure scipy is installed.") from e

        # 4) Cm
        cm_full = estimate_Cm_full(alpha_full, df)

        # 5) Save CSV
        out_df = pd.DataFrame({
            'alpha_deg': alpha_full,
            'CL': CL_full,
            'CD': CD_full,
            'CM': cm_full
        })
        out_df.to_csv(out_csv, index=False)
        if self.verbose:
            print(f"[save] Wrote 0-360 polar to '{out_csv}'")

        # 6) Plot
        fig, axs = plt.subplots(3, 1, figsize=(8, 9), constrained_layout=True)
        axs[0].plot(alpha_full, CL_full)
        axs[0].set_ylabel("CL")
        axs[0].set_xlim(0, 360)
        axs[0].grid(True)

        axs[1].plot(alpha_full, CD_full)
        axs[1].set_ylabel("CD")
        axs[1].set_xlim(0, 360)
        axs[1].grid(True)

        axs[2].plot(alpha_full, cm_full)
        axs[2].set_ylabel("CM")
        axs[2].set_xlim(0, 360)
        axs[2].set_xlabel("alpha (deg)")
        axs[2].grid(True)

        plt.suptitle(f"0-360° Polar: {airfoil}")
        plt.savefig(plot_fig, dpi=180)
        if self.verbose:
            print(f"[plot] Saved plot to '{plot_fig}'")
        plt.show()

        return out_df

# -------------------------
# Example script usage
# -------------------------
def main():
    # Example parameters — change to suit your case
    xfoil_exec = xfoil_path  # or "C:\\xfoil\\xfoil.exe" on Windows
    airfoil = "NACA2412"  # or "mh114.dat"
    alphas = np.linspace(-8, 15, 24)  # XFOIL sweep range (pre-stall region)
    Re = 5e5
    AR = 6.0

    x360 = XfoilPolar360(xfoil_path=xfoil_exec, verbose=True)
    out = x360.make_polar_and_extrapolate(airfoil=airfoil,
                                          alphas=alphas,
                                          Re=Re,
                                          Mach=0.0,
                                          AR=AR,
                                          polar_filename="polar.txt",
                                          out_csv="polar_0_360.csv",
                                          plot_fig="polar_0_360.png")
    print("Done. Saved CSV head:")
    print(out.head())

if __name__ == "__main__":
    main()
