import math
import numpy as np
from scipy.interpolate import interp1d

# ============ Load Airfoil Data from CSV Files ============
def load_airfoil_data(cl_csv_path, cd_csv_path):
    """Load Cl and Cd data from CSV files and return interpolation functions."""
    try:
        # Load Cl data - handle both comma and semicolon delimiters, skip header rows
        cl_data = np.genfromtxt(cl_csv_path, delimiter=';', skip_header=3)
        alpha_cl_deg = cl_data[:, 0]  # Keep in degrees for now
        alpha_cl = np.radians(alpha_cl_deg)  # Convert to radians
        Cl_values = cl_data[:, 1]
        
        # Load Cd data - handle both comma and semicolon delimiters, skip header rows
        cd_data = np.genfromtxt(cd_csv_path, delimiter=';', skip_header=3)
        alpha_cd_deg = cd_data[:, 0]  # Keep in degrees for now
        alpha_cd = np.radians(alpha_cd_deg)  # Convert to radians
        Cd_values = cd_data[:, 1]
        
        # Create interpolation functions with linear interpolation
        # Using bounds_error=False allows extrapolation at boundaries
        cl_func = interp1d(alpha_cl, Cl_values, kind='linear', bounds_error=False, 
                          fill_value=(Cl_values[0], Cl_values[-1]))
        cd_func = interp1d(alpha_cd, Cd_values, kind='linear', bounds_error=False,
                          fill_value=(Cd_values[0], Cd_values[-1]))
        
        print(f"✓ Successfully loaded 360° polar data from CSV files")
        print(f"  Cl: {len(Cl_values)} points, α = {alpha_cl_deg[0]:.1f}° to {alpha_cl_deg[-1]:.1f}°")
        print(f"      Range: {Cl_values.min():.4f} to {Cl_values.max():.4f}")
        print(f"  Cd: {len(Cd_values)} points, α = {alpha_cd_deg[0]:.1f}° to {alpha_cd_deg[-1]:.1f}°")
        print(f"      Range: {Cd_values.min():.6f} to {Cd_values.max():.6f}")
        
        return cl_func, cd_func
    except Exception as e:
        print(f"⚠ Error loading airfoil data: {e}")
        print(f"  Falling back to theoretical models...")
        import traceback
        traceback.print_exc()
        return None, None

# Default paths for CSV files (relative to script location)
# DEFAULT_CL_CSV = 'naca2412_lift_coeff.csv'
# DEFAULT_CD_CSV = 'naca2412_drag_coeff.csv'
DEFAULT_CL_CSV = 'naca2412_360_polar_lift.csv'
DEFAULT_CD_CSV = 'naca2412_360_polar_drag.csv'

# Global interpolation functions (will be set by load_airfoil_data)
_cl_interp = None
_cd_interp = None

# ============ Prandtl's Loss Factors ============
def prandtl_tip_loss(r, R, B, phi):
    """Prandtl's tip loss factor."""
    f = 0.5 * B * (1.0 - r/R) / (r/R * abs(math.sin(phi)))
    return (2.0/math.pi) * math.acos(min(1.0, math.exp(-f)))

def prandtl_root_loss(r, R_hub, B, phi):
    """Prandtl's root loss factor."""
    g = 0.5 * B * (r/R_hub - 1.0) / (r/R_hub * abs(math.sin(phi)))
    return (2.0/math.pi) * math.acos(min(1.0, math.exp(-g)))

# ============ Aerodynamic models with CSV support ============
def cl_from_alpha(alpha, cl_func=None, Cl_alpha=2*math.pi):
    """
    Get lift coefficient from angle of attack.
    Uses CSV interpolation if available, otherwise uses linear model.
    
    Parameters:
    -----------
    alpha : float
        Angle of attack (radians)
    cl_func : callable
        Interpolation function from CSV data
    Cl_alpha : float
        Lift slope for linear model (default = 2π per rad)
    """
    global _cl_interp
    
    # Use provided function or global interpolation function
    func = cl_func if cl_func is not None else _cl_interp
    
    if func is not None:
        # Use CSV data via interpolation
        return float(func(alpha))
    else:
        # Fall back to linear model
        return Cl_alpha * alpha

def cd_from_alpha(alpha, cd_func=None, Cd0=0.008, k=0.02):
    """
    Get drag coefficient from angle of attack.
    Uses CSV interpolation if available, otherwise uses quadratic model.
    
    Parameters:
    -----------
    alpha : float
        Angle of attack (radians)
    cd_func : callable
        Interpolation function from CSV data
    Cd0 : float
        Parasite drag coefficient (for fallback model)
    k : float
        Drag due to lift parameter (for fallback model)
    """
    global _cd_interp
    
    # Use provided function or global interpolation function
    func = cd_func if cd_func is not None else _cd_interp
    
    if func is not None:
        # Use CSV data via interpolation
        return float(func(alpha))
    else:
        # Fall back to quadratic model
        return Cd0 + k * (alpha**2)

# ---------------- Helper functions ----------------
def compute_gamma(Cl, Cd):
    if abs(Cl) < 1e-12:
        return math.atan2(Cd, 1e-12)
    return math.atan2(Cd, Cl)

def F_func(phi, sigma, Cl, Cd, F_tip, F_root):
    """Modified F function with tip and root loss factors."""
    sec = 1.0 / math.cos(phi)
    gamma = compute_gamma(Cl, Cd)
    F = F_tip * F_root  # Combined loss factor
    sin_phi = math.sin(phi)
    cosec_phi = 1.0 / sin_phi
    cos_gamma = math.cos(gamma)
    sec_gamma = 1.0 / cos_gamma
    cos_gamma_phi = math.cos(gamma + phi)
    f = sin_phi - 0.25 * sigma * sec_gamma * Cl * cosec_phi * cos_gamma_phi
    return f

def G_func(phi, sigma, Cl, Cd, F_tip, F_root):
    """Modified G function with tip and root loss factors."""
    sin_phi = math.sin(phi)
    cos_phi = math.cos(phi)
    gamma = compute_gamma(Cl, Cd)
    cos_gamma = math.cos(gamma)
    sec_gamma = 1.0 / cos_gamma
    csc_phi = 1.0 / sin_phi  # CRITICAL: cosecant term
    sin_gamma_phi = math.sin(gamma + phi)
    F = F_tip * F_root  # Combined loss factor
    g = cos_phi + 0.25 * sigma * sec_gamma * Cl * csc_phi * sin_gamma_phi
    return g

# ---------------- Glauert Correction ----------------
def glauert_correction(a):
    """Glauert empirical correction for heavily loaded rotors."""
    a_c = 0.2  # Critical induction factor
    if a <= a_c:
        return a
    return 0.5 * (2 + math.sqrt(4 - 4*(1 - 2*a_c)))
    # return 0.25 * (5 - (3 * a))
    # return (a_c/a)*(2 - (a_c/a))

# ---------------- φ Solver (Modified Regula Falsi) ----------------
def solve_phi_regula_falsi(r, V, Omega, sigma, theta_local, Cl_alpha_func, Cd_func,
                          R, R_hub, B, alpha_zero_lift=0.0, pitch_beta=0.0,
                          max_iter=80, tol=1e-8):
    """Find inflow angle φ using a bracketed root method with tip/root loss."""
    left = 1e-6
    right = math.pi/2 - 1e-6

    def g_of_phi(phi):
        # Calculate loss factors
        F_tip = prandtl_tip_loss(r, R, B, phi)
        F_root = prandtl_root_loss(r, R_hub, B, phi)
        
        alpha = (theta_local + pitch_beta - alpha_zero_lift) - phi
        Cl = Cl_alpha_func(alpha)
        Cd = Cd_func(alpha)
        F = F_func(phi, sigma, Cl, Cd, F_tip, F_root)
        G = G_func(phi, sigma, Cl, Cd, F_tip, F_root)
        
        denom = r * Omega if abs(r * Omega) > 1e-12 else 1e-12
        return math.sin(phi) * (F - (V / denom) * G)

    f_left = g_of_phi(left)
    f_right = g_of_phi(right)

    # Fallback if bracket not valid
    if f_left * f_right > 0:
        for _ in range(max_iter):
            mid = 0.5 * (left + right)
            f_mid = g_of_phi(mid)
            if abs(f_mid) < tol:
                return mid
            if f_left * f_mid <= 0:
                right, f_right = mid, f_mid
            else:
                left, f_left = mid, f_mid
        return 0.5 * (left + right)

    phi = None
    for iteration in range(max_iter):
        denom = f_right - f_left
        if abs(denom) < 1e-15:  # Prevent division by zero
            return 0.5 * (left + right)
        
        phi_new = (left * f_right - right * f_left) / denom
        
        # Ensure phi_new is within bounds
        if phi_new <= left or phi_new >= right:
            phi_new = 0.5 * (left + right)
        
        f_new = g_of_phi(phi_new)
        
        if abs(f_new) < tol:
            return phi_new
        
        if f_left * f_new < 0:
            right, f_right = phi_new, f_new
        else:
            left, f_left = phi_new, f_new
        
        phi = phi_new
        
        # Check for convergence
        if abs(right - left) < tol:
            return phi
    
    return phi if phi is not None else 0.5 * (left + right)

# ============ Elemental Forces ============
def elemental_forces(r, c, B, rho, V, Omega, theta_local, pitch_beta, alpha0,
                    Cl_alpha=2*math.pi, Cd0=0.008, Cd_k=0.02, R_hub=0.05, R=0.5,
                    cl_func=None, cd_func=None):
    """
    Compute elemental thrust & torque with tip/root loss effects.
    
    Parameters:
    -----------
    r, c, B, rho, V, Omega, theta_local, pitch_beta, alpha0 : float
        BEMT input parameters
    Cl_alpha, Cd0, Cd_k : float
        Fallback model parameters (used if CSV data not available)
    R_hub, R : float
        Hub and tip radius
    cl_func, cd_func : callable
        Interpolation functions from CSV data (optional)
    """
    global _cl_interp, _cd_interp
    
    sigma = B * c / (2.0 * math.pi * r)
    
    # Use provided functions or global interpolation functions
    cl_fn = lambda alpha: cl_from_alpha(alpha, cl_func=cl_func or _cl_interp, Cl_alpha=Cl_alpha)
    cd_fn = lambda alpha: cd_from_alpha(alpha, cd_func=cd_func or _cd_interp, Cd0=Cd0, k=Cd_k)

    phi = solve_phi_regula_falsi(r, V, Omega, sigma, theta_local, cl_fn, cd_fn,
                                R, R_hub, B, alpha_zero_lift=alpha0, pitch_beta=pitch_beta)

    # Calculate loss factors
    F_tip = prandtl_tip_loss(r, R, B, phi)
    F_root = prandtl_root_loss(r, R_hub, B, phi)
    
    alpha = (theta_local + pitch_beta - alpha0) - phi
    Cl = cl_fn(alpha)
    Cd = cd_fn(alpha)
   

    # Calculate induction factors with Glauert correction
    F = F_func(phi, sigma, Cl, Cd, F_tip, F_root)
    G = G_func(phi, sigma, Cl, Cd, F_tip, F_root)
    a = (1.0 - (math.sin(phi) / F)) if abs(F) > 1e-12 else 0.0
    # a = glauert_correction(a)  # Apply Glauert correction
    b = (math.cos(phi) / G - 1.0) if abs(G) > 1e-12 else 0.0
    # TRY: Use V(1-a) convention instead of V(1+a)
    V_axial = V * (1 - a)
    V_s = 2*V_axial - V
    V_tangential = r * Omega * (1 + b)

    Vr = math.sqrt(V_axial**2 + (V_tangential)**2)

    # Apply loss factors to thrust and torque
    F_total = F_tip * F_root
    # FIX: Drag contributes positively to thrust (not negatively)
    dT = 0.5 * B * c * rho * (Vr**2) * (Cl * math.cos(phi) + Cd * math.sin(phi)) * F_total
    dQ = 0.5 * B * c * rho * (Vr**2) * (Cl * math.sin(phi) + Cd * math.cos(phi)) * r * F_total

    diagnostics = dict(
        r=r, phi=phi, alpha=alpha, Cl=Cl, Cd=Cd, sigma=sigma,
        a=a, b=b, F_tip=F_tip, F_root=F_root, dT=dT, dQ=dQ, V_s=V_s
    )
    return dT, dQ, diagnostics

# ============ Integration Across Blade ============
def integrate_blade(B, chord_fun, theta_fun, R, R_hub, N, rho, V, rpm,
                   pitch_beta=0.0, alpha0_fun=None,
                   Cl_alpha=2*math.pi, Cd0=0.008, Cd_k=0.02,
                   cl_func=None, cd_func=None):
    """
    Integrate forces along blade with improved numerical integration.
    
    Parameters:
    -----------
    cl_func, cd_func : callable
        Interpolation functions from CSV data (optional)
    """
    global _cl_interp, _cd_interp
    
    Omega = 2.0 * math.pi * rpm / 60.0
    
    # Cosine spacing for better resolution near tips
    beta = np.linspace(0, math.pi, N+1)
    nodes = R_hub + (R - R_hub) * (1 - np.cos(beta))/2
    
    dT_nodes, dQ_nodes, diags = [], [], []

    for r in nodes:
        if r <= R_hub + 1e-12:
            dT_nodes.append(0.0)
            dQ_nodes.append(0.0)
            diags.append({})
            continue
            
        c = chord_fun(r)
        theta_local = theta_fun(r)
        alpha0 = alpha0_fun(r) if alpha0_fun is not None else 0.0
        
        dT, dQ, diagnostics = elemental_forces(
            r, c, B, rho, V, Omega, theta_local, pitch_beta, alpha0,
            Cl_alpha=Cl_alpha, Cd0=Cd0, Cd_k=Cd_k, R_hub=R_hub, R=R,
            cl_func=cl_func or _cl_interp, cd_func=cd_func or _cd_interp
        )
        
        dT_nodes.append(dT)
        dQ_nodes.append(dQ)
        diags.append(diagnostics)

    # Trapezoidal integration
    T_total, Q_total = 0.0, 0.0
    for i in range(N):
        dr = nodes[i+1] - nodes[i]
        T_total += 0.5 * (dT_nodes[i] + dT_nodes[i+1]) * dr
        Q_total += 0.5 * (dQ_nodes[i] + dQ_nodes[i+1]) * dr

    return T_total, Q_total, diags, nodes

# ============ Example Run ============
if __name__ == "__main__":
    import os
    
    # Try to load CSV files
    script_dir = os.path.dirname(os.path.abspath(__file__))
    cl_csv = os.path.join(script_dir, DEFAULT_CL_CSV)
    cd_csv = os.path.join(script_dir, DEFAULT_CD_CSV)
    
    print("="*60)
    print("BEMT ANALYSIS WITH CSV AIRFOIL DATA")
    print("="*60)
    
    # Load airfoil data
    _cl_interp, _cd_interp = load_airfoil_data(cl_csv, cd_csv)
    
    # Propeller parameters
    B = 3
    R = 0.455
    R_hub = 0.05
    rho = 1.225
    rpm = 600        # Nominal RPM
    V = 5.0           # Low velocity for propeller mode
    N = 20
    A = math.pi * (R**2 )

    # Blade geometry functions
    chord_fun = lambda r: 0.07
    # theta_fun = lambda r: math.radians(20.0) * (1 - (r - R_hub)/(R - R_hub))  # 20° twist
    theta_fun = lambda r: 0.0
    alpha0_fun = lambda r: 0.0

    print(f"\nPropeller Configuration:")
    print(f"  Blades: {B}")
    print(f"  Radius: {R} m, Hub: {R_hub} m")
    print(f"  Chord: 0.2 m (constant)")
    print(f"  Twist: Linear {math.degrees(theta_fun(R_hub)):.1f}° (root) to {math.degrees(theta_fun(R)):.1f}° (tip)")
    
    print(f"\nOperating Condition:")
    print(f"  Velocity: {V} m/s")
    print(f"  RPM: {rpm}")
    print(f"  Blade elements: {N}")
    
    # Run analysis with CSV data
    T, Q, diags, nodes = integrate_blade(
        B, chord_fun, theta_fun, R, R_hub, N, rho, V, rpm,
        pitch_beta=math.radians(0.0),
        alpha0_fun=alpha0_fun,
        Cl_alpha=2*math.pi, Cd0=0.008, Cd_k=0.02,
        cl_func=_cl_interp, cd_func=_cd_interp
    )

    print(f"\n" + "="*60)
    print("RESULTS:")
    print("="*60)
    print(f"Total Thrust = {T:.4f} N")
    print(f"Total Torque = {Q:.4f} N·m")
    Omega = 2*math.pi*rpm/60
    # w = math.sqrt(2*T/(rho*A))
    P_in = Q*Omega
    P_out = T*V
    eta = 100*P_out/P_in if P_in > 0 else 0
    print(f"Input Power = {P_in:.2f} W")
    print(f"Output Power = {P_out:.2f} W")
    print(f"Efficiency = {eta:.1f}%")
    # print(f"Wake induced velocity = {w:.2f} m/s")

    # Print mid-span diagnostics
    mid = len(nodes)//2
    print(f"\nMid-span station diagnostics (r/R = {nodes[mid]/R:.3f}):")
    print("-"*60)
    for key, value in diags[mid].items():
        if isinstance(value, float):
            print(f"{key:>10s}: {value:>12.6f}")
        else:
            print(f"{key}: {value}")
    
    # Test multiple velocities to demonstrate CSV integration
    print(f"\n" + "="*60)
    print("MULTI-VELOCITY TEST (CSV DATA INTEGRATION):")
    print("="*60)
    test_velocities = [5.0, 8.0, 10.0, 12.0, 15.0]
    for V_test in test_velocities:
        T_test, Q_test, _, _ = integrate_blade(
            B, chord_fun, theta_fun, R, R_hub, N, rho, V_test, rpm,
            pitch_beta=math.radians(0.0), alpha0_fun=alpha0_fun,
            Cl_alpha=2*math.pi, Cd0=0.008, Cd_k=0.02,
            cl_func=_cl_interp, cd_func=_cd_interp
        )
        P_in_test = Q_test * Omega
        P_out_test = T_test * V_test
        eta_test = 100*P_out_test/P_in_test if P_in_test > 0 else 0
        print(f"V={V_test:5.1f} m/s: T={T_test:7.3f} N, Q={Q_test:9.3f} N·m, η={eta_test:6.1f}%")