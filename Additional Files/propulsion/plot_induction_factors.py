import numpy as np
import matplotlib.pyplot as plt
from bemt_rmit_improved import integrate_blade
import math

def plot_induction_factors_vs_radial_distance(B, chord_fun, theta_fun, R, R_hub, N, 
                                             rho, V, rpm, pitch_beta=0.0, 
                                             alpha0_fun=None, Cl_alpha=2*math.pi, 
                                             Cd0=0.008, Cd_k=0.02):
    """
    Plot induction factors (a and a') and inflow angle (φ) vs radial distance.
    
    Parameters:
    -----------
    B : int
        Number of blades
    chord_fun : callable
        Function that returns chord as f(r)
    theta_fun : callable
        Function that returns blade pitch angle as f(r)
    R : float
        Blade radius (m)
    R_hub : float
        Hub radius (m)
    N : int
        Number of blade elements
    rho : float
        Air density (kg/m³)
    V : float
        Freestream velocity (m/s)
    rpm : float
        Rotational speed (rpm)
    pitch_beta : float
        Collective pitch angle (rad)
    alpha0_fun : callable
        Function for angle of zero lift
    Cl_alpha : float
        Lift curve slope
    Cd0 : float
        Parasite drag coefficient
    Cd_k : float
        Drag due to lift parameter
    """
    
    # Run BEMT calculation
    T, Q, diags, nodes = integrate_blade(
        B, chord_fun, theta_fun, R, R_hub, N, rho, V, rpm,
        pitch_beta=pitch_beta, alpha0_fun=alpha0_fun,
        Cl_alpha=Cl_alpha, Cd0=Cd0, Cd_k=Cd_k
    )
    
    # Extract data from diagnostics
    r_values = np.array([d['r'] for d in diags if d])  # Absolute radius
    r_R = r_values / R  # Non-dimensional radius
    
    a = np.array([d['a'] for d in diags if d])  # Axial induction factor
    b = np.array([d['b'] for d in diags if d])  # Tangential induction factor
    phi = np.array([d['phi'] for d in diags if d])  # Inflow angle (radians)
    
    # Create figure with subplots
    fig = plt.figure(figsize=(14, 10))
    
    # Plot 1: Induction Factors vs r/R
    ax1 = plt.subplot(2, 2, 1)
    ax1.plot(r_R, a, 'b-', linewidth=2.5, marker='o', markersize=5, label='Axial induction (a)')
    ax1.plot(r_R, b, 'r--', linewidth=2.5, marker='s', markersize=5, label='Tangential induction (a\')')
    ax1.set_xlabel('Radial Distance (r/R)', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Induction Factor', fontsize=11, fontweight='bold')
    ax1.set_title('Induction Factors vs Radial Distance', fontsize=12, fontweight='bold')
    ax1.legend(fontsize=10, loc='best')
    ax1.grid(True, alpha=0.3, linestyle='--')
    ax1.set_xlim([0, 1])
    
    # Plot 2: Inflow Angle vs r/R (in degrees)
    ax2 = plt.subplot(2, 2, 2)
    phi_deg = np.degrees(phi)
    ax2.plot(r_R, phi_deg, 'g-', linewidth=2.5, marker='^', markersize=5)
    ax2.set_xlabel('Radial Distance (r/R)', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Inflow Angle (degrees)', fontsize=11, fontweight='bold')
    ax2.set_title('Inflow Angle (φ) vs Radial Distance', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3, linestyle='--')
    ax2.set_xlim([0, 1])
    
    # Plot 3: Combined plot of all three factors (normalized)
    ax3 = plt.subplot(2, 2, 3)
    ax3.plot(r_R, a, 'b-', linewidth=2.5, marker='o', markersize=5, label='a (axial)')
    ax3.plot(r_R, b, 'r--', linewidth=2.5, marker='s', markersize=5, label='a\' (tangential)')
    ax3_twin = ax3.twinx()  # Create secondary y-axis
    ax3_twin.plot(r_R, phi_deg, 'g-', linewidth=2.5, marker='^', markersize=5, label='φ (degrees)')
    
    ax3.set_xlabel('Radial Distance (r/R)', fontsize=11, fontweight='bold')
    ax3.set_ylabel('Induction Factors', fontsize=11, fontweight='bold', color='black')
    ax3_twin.set_ylabel('Inflow Angle (degrees)', fontsize=11, fontweight='bold', color='green')
    ax3.set_title('All Parameters vs Radial Distance', fontsize=12, fontweight='bold')
    ax3.tick_params(axis='y', labelcolor='black')
    ax3_twin.tick_params(axis='y', labelcolor='green')
    ax3.grid(True, alpha=0.3, linestyle='--')
    ax3.set_xlim([0, 1])
    
    # Add legends from both axes
    lines1, labels1 = ax3.get_legend_handles_labels()
    lines2, labels2 = ax3_twin.get_legend_handles_labels()
    ax3.legend(lines1 + lines2, labels1 + labels2, fontsize=10, loc='upper right')
    
    # Plot 4: Statistics table
    ax4 = plt.subplot(2, 2, 4)
    ax4.axis('off')
    
    # Calculate statistics
    Omega = 2 * math.pi * rpm / 60
    V_tip = Omega * R
    J = V / (rpm / 60 * 2 * R)  # Advance ratio
    
    # Find values at different radii
    idx_root = 1 if len(r_R) > 1 else 0
    idx_mid = len(r_R) // 2
    idx_tip = -1
    
    # Create statistics text
    stats_text = f"""
BEMT CALCULATION SUMMARY
{'='*40}

Operating Conditions:
  • Freestream Velocity (V):     {V:.2f} m/s
  • RPM:                          {rpm:.0f}
  • Collective Pitch (β):          {np.degrees(pitch_beta):.2f}°
  • Advance Ratio (J):            {J:.3f}

Blade Geometry:
  • Number of Blades (B):         {B}
  • Blade Radius (R):             {R:.3f} m
  • Hub Radius (R_hub):           {R_hub:.3f} m
  • Root Chord:                   {chord_fun(R_hub):.4f} m
  • Tip Chord:                    {chord_fun(R):.4f} m

Performance:
  • Total Thrust:                 {T:.3f} N
  • Total Torque:                 {Q:.3f} N·m
  • Input Power:                  {Q*Omega:.2f} W
  • Output Power:                 {T*V:.2f} W
  • Efficiency:                   {100*T*V/(Q*Omega):.1f}%

Induction Factors (Root/Mid/Tip):
  • Axial (a):            {a[idx_root]:.4f} / {a[idx_mid]:.4f} / {a[idx_tip]:.4f}
  • Tangential (a'):      {b[idx_root]:.4f} / {b[idx_mid]:.4f} / {b[idx_tip]:.4f}

Inflow Angle (Root/Mid/Tip):
  • φ (degrees):          {phi_deg[idx_root]:.2f}° / {phi_deg[idx_mid]:.2f}° / {phi_deg[idx_tip]:.2f}°
  • φ (radians):          {phi[idx_root]:.4f} / {phi[idx_mid]:.4f} / {phi[idx_tip]:.4f}
"""
    
    ax4.text(0.05, 0.95, stats_text, transform=ax4.transAxes, 
             fontsize=9, verticalalignment='top', family='monospace',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    return fig

if __name__ == "__main__":
    # Example propeller parameters
    B = 3  # Number of blades
    R = 0.5  # Blade radius (m)
    R_hub = 0.05  # Hub radius (m)
    rho = 1.225  # Air density (kg/m³)
    N = 40  # Number of blade elements
    
    # Blade geometry functions
    # Constant chord
    chord_fun = lambda r: 0.04  # 4 cm chord
    
    # Linear twist (higher pitch at root, lower at tip)
    theta_fun = lambda r: math.radians(20.0) * (1 - (r - R_hub)/(R - R_hub))
    
    # Operating condition 1: Hover-like condition
    print("Generating plot for Condition 1: V=5 m/s, RPM=4000")
    fig1 = plot_induction_factors_vs_radial_distance(
        B, chord_fun, theta_fun, R, R_hub, N, rho,
        V=5.0, rpm=4000, pitch_beta=math.radians(0)
    )
    fig1.savefig('induction_factors_condition1.png', dpi=300, bbox_inches='tight')
    
    # Operating condition 2: Forward flight condition
    print("Generating plot for Condition 2: V=10 m/s, RPM=4000")
    fig2 = plot_induction_factors_vs_radial_distance(
        B, chord_fun, theta_fun, R, R_hub, N, rho,
        V=10.0, rpm=4000, pitch_beta=math.radians(0)
    )
    fig2.savefig('induction_factors_condition2.png', dpi=300, bbox_inches='tight')
    
    # Operating condition 3: High speed condition
    print("Generating plot for Condition 3: V=15 m/s, RPM=4000")
    fig3 = plot_induction_factors_vs_radial_distance(
        B, chord_fun, theta_fun, R, R_hub, N, rho,
        V=15.0, rpm=4000, pitch_beta=math.radians(0)
    )
    fig3.savefig('induction_factors_condition3.png', dpi=300, bbox_inches='tight')
    
    plt.show()
    
    print("\nPlots saved as PNG files:")
    print("  - induction_factors_condition1.png")
    print("  - induction_factors_condition2.png")
    print("  - induction_factors_condition3.png")
