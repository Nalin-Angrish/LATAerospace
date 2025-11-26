"""
Interactive inline plotting of induction factors and inflow angles
Shows detailed annotations and hover information
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from bemt_rmit_improved import integrate_blade
import math

def create_annotated_single_plot(B, chord_fun, theta_fun, R, R_hub, N, rho, V, rpm):
    """
    Create a detailed annotated plot with enhanced visualization
    """
    # Run BEMT
    T, Q, diags, nodes = integrate_blade(
        B, chord_fun, theta_fun, R, R_hub, N, rho, V, rpm
    )
    
    # Extract data
    r_values = np.array([d['r'] for d in diags if d])
    r_R = r_values / R
    a = np.array([d['a'] for d in diags if d])
    b = np.array([d['b'] for d in diags if d])
    phi = np.array([d['phi'] for d in diags if d])
    alpha = np.array([d['alpha'] for d in diags if d])
    Cl = np.array([d['Cl'] for d in diags if d])
    Cd = np.array([d['Cd'] for d in diags if d])
    
    phi_deg = np.degrees(phi)
    alpha_deg = np.degrees(alpha)
    
    # Create figure with detailed layout
    fig = plt.figure(figsize=(18, 12))
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
    
    # ===== ROW 1: Induction Factors =====
    
    # Plot 1: Axial induction
    ax1 = fig.add_subplot(gs[0, 0])
    line1a = ax1.plot(r_R, a, 'b-', linewidth=3, label='a (axial)')
    ax1.fill_between(r_R, a, alpha=0.15, color='blue')
    ax1.scatter(r_R[::5], a[::5], color='blue', s=80, zorder=5)
    ax1.set_ylabel('Axial Induction (a)', fontsize=11, fontweight='bold')
    ax1.set_title('Axial Induction Factor', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3, linestyle='--')
    ax1.set_xlim([0, 1])
    ax1.text(0.5, 0.95, f'max={a.max():.4f}, min={a.min():.4f}', 
             transform=ax1.transAxes, ha='center', va='top',
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.7), fontsize=9)
    
    # Plot 2: Tangential induction
    ax2 = fig.add_subplot(gs[0, 1])
    line1b = ax2.plot(r_R, b, 'r-', linewidth=3, label='a\' (tangential)')
    ax2.fill_between(r_R, b, alpha=0.15, color='red')
    ax2.scatter(r_R[::5], b[::5], color='red', s=80, zorder=5)
    ax2.set_ylabel('Tangential Induction (a\')', fontsize=11, fontweight='bold')
    ax2.set_title('Tangential Induction Factor', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3, linestyle='--')
    ax2.set_xlim([0, 1])
    ax2.text(0.5, 0.95, f'max={b.max():.4f}, min={b.min():.4f}', 
             transform=ax2.transAxes, ha='center', va='top',
             bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.7), fontsize=9)
    
    # Plot 3: Combined induction
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.plot(r_R, a, 'b-', linewidth=2.5, label='a', marker='o', markersize=4)
    ax3.plot(r_R, b, 'r-', linewidth=2.5, label='a\'', marker='s', markersize=4)
    ax3.set_ylabel('Induction Factor', fontsize=11, fontweight='bold')
    ax3.set_title('Combined Induction Factors', fontsize=12, fontweight='bold')
    ax3.legend(fontsize=10, loc='best')
    ax3.grid(True, alpha=0.3, linestyle='--')
    ax3.set_xlim([0, 1])
    
    # ===== ROW 2: Angles =====
    
    # Plot 4: Inflow angle
    ax4 = fig.add_subplot(gs[1, 0])
    ax4.plot(r_R, phi_deg, 'g-', linewidth=3)
    ax4.fill_between(r_R, phi_deg, alpha=0.15, color='green')
    ax4.scatter(r_R[::5], phi_deg[::5], color='green', s=80, zorder=5)
    ax4.set_ylabel('Inflow Angle (degrees)', fontsize=11, fontweight='bold')
    ax4.set_title('Inflow Angle φ', fontsize=12, fontweight='bold')
    ax4.grid(True, alpha=0.3, linestyle='--')
    ax4.set_xlim([0, 1])
    ax4.text(0.5, 0.95, f'max={phi_deg.max():.2f}°, min={phi_deg.min():.2f}°', 
             transform=ax4.transAxes, ha='center', va='top',
             bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7), fontsize=9)
    
    # Plot 5: Angle of attack
    ax5 = fig.add_subplot(gs[1, 1])
    ax5.plot(r_R, alpha_deg, 'm-', linewidth=3)
    ax5.fill_between(r_R, alpha_deg, alpha=0.15, color='magenta')
    ax5.scatter(r_R[::5], alpha_deg[::5], color='magenta', s=80, zorder=5)
    ax5.set_ylabel('Angle of Attack (degrees)', fontsize=11, fontweight='bold')
    ax5.set_title('Angle of Attack α', fontsize=12, fontweight='bold')
    ax5.grid(True, alpha=0.3, linestyle='--')
    ax5.set_xlim([0, 1])
    ax5.axhline(y=0, color='k', linestyle='--', alpha=0.3)
    ax5.text(0.5, 0.95, f'max={alpha_deg.max():.2f}°, min={alpha_deg.min():.2f}°', 
             transform=ax5.transAxes, ha='center', va='top',
             bbox=dict(boxstyle='round', facecolor='plum', alpha=0.7), fontsize=9)
    
    # Plot 6: Aerodynamic coefficients
    ax6 = fig.add_subplot(gs[1, 2])
    ax6_twin = ax6.twinx()
    line6a = ax6.plot(r_R, Cl, 'b-', linewidth=2.5, label='Cl', marker='o', markersize=4)
    line6b = ax6_twin.plot(r_R, Cd, 'r-', linewidth=2.5, label='Cd', marker='s', markersize=4)
    ax6.set_ylabel('Lift Coefficient (Cl)', fontsize=11, fontweight='bold', color='blue')
    ax6_twin.set_ylabel('Drag Coefficient (Cd)', fontsize=11, fontweight='bold', color='red')
    ax6.set_title('Aerodynamic Coefficients', fontsize=12, fontweight='bold')
    ax6.tick_params(axis='y', labelcolor='blue')
    ax6_twin.tick_params(axis='y', labelcolor='red')
    ax6.grid(True, alpha=0.3, linestyle='--')
    ax6.set_xlim([0, 1])
    
    # ===== ROW 3: Summary Information =====
    
    # Plot 7: Velocity diagram
    ax7 = fig.add_subplot(gs[2, 0])
    # Show velocity components at mid-span
    idx_mid = len(r_R) // 2
    r_mid = r_values[idx_mid]
    phi_mid = phi[idx_mid]
    a_mid = a[idx_mid]
    b_mid = b[idx_mid]
    
    Omega = 2 * math.pi * rpm / 60
    V_axial = V * (1 - a_mid)
    V_tangential = Omega * r_mid * (1 + b_mid)
    
    # Vector plot
    ax7.quiver(0, 0, V, 0, angles='xy', scale_units='xy', scale=1, color='blue', 
              width=0.01, label='Freestream V')
    ax7.quiver(0, 0, V_axial, V_tangential, angles='xy', scale_units='xy', scale=1, 
              color='red', width=0.01, label='Relative velocity')
    ax7.quiver(0, 0, 0, Omega*r_mid, angles='xy', scale_units='xy', scale=1, 
              color='green', width=0.01, label='Blade velocity')
    
    ax7.set_xlim([-2, V+2])
    ax7.set_ylim([-2, max(Omega*r_mid, V_tangential)+2])
    ax7.set_aspect('equal')
    ax7.grid(True, alpha=0.3)
    ax7.set_xlabel('Axial Velocity (m/s)', fontsize=10, fontweight='bold')
    ax7.set_ylabel('Tangential Velocity (m/s)', fontsize=10, fontweight='bold')
    ax7.set_title(f'Velocity Diagram (r/R={r_R[idx_mid]:.2f})', fontsize=11, fontweight='bold')
    ax7.legend(fontsize=9, loc='upper left')
    
    # Plot 8: Performance metrics
    ax8 = fig.add_subplot(gs[2, 1])
    ax8.axis('off')
    
    Omega = 2 * math.pi * rpm / 60
    J = V / (rpm/60 * 2*R)
    P_in = Q * Omega
    P_out = T * V
    eta = 100 * P_out / P_in if P_in > 0 else 0
    
    metrics_text = f"""
PERFORMANCE SUMMARY
{'─'*30}
Velocity: {V:.2f} m/s
RPM: {rpm:.0f}
Adv. Ratio (J): {J:.4f}

Thrust: {T:.4f} N
Torque: {Q:.4f} N·m
Input Power: {P_in:.2f} W
Output Power: {P_out:.2f} W
Efficiency: {eta:.1f}%

Blade: R={R:.3f}m, Nhub={R_hub:.3f}m
N_blades={B}, Elements={N}

Avg Induction:
  a_avg: {a.mean():.4f}
  a'_avg: {b.mean():.4f}
"""
    ax8.text(0.1, 0.95, metrics_text, transform=ax8.transAxes,
            fontsize=9, verticalalignment='top', family='monospace',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    
    # Plot 9: Data table at key radii
    ax9 = fig.add_subplot(gs[2, 2])
    ax9.axis('off')
    
    # Create table
    table_data = []
    indices = [0, len(r_R)//4, len(r_R)//2, 3*len(r_R)//4, -1]
    
    table_text = f"""
RADIAL STATION DATA
{'─'*45}
r/R    a       a'      φ°    α°    Cl
{'─'*45}
"""
    for i in indices:
        if i < len(r_R):
            table_text += f"{r_R[i]:.3f}  {a[i]:.4f}  {b[i]:.4f}  "
            table_text += f"{phi_deg[i]:6.2f}  {alpha_deg[i]:6.2f}  {Cl[i]:.4f}\n"
    
    ax9.text(0.05, 0.95, table_text, transform=ax9.transAxes,
            fontsize=8, verticalalignment='top', family='monospace',
            bbox=dict(boxstyle='round', facecolor='lightcyan', alpha=0.8))
    
    # Add overall title
    fig.suptitle(f'BEMT Induction Analysis: V={V} m/s, RPM={rpm}, 3-Blade Propeller',
                fontsize=14, fontweight='bold', y=0.995)
    
    # Add x-label to bottom row
    for ax in [ax7, ax8, ax9]:
        if hasattr(ax, 'set_xlabel'):
            pass  # Skip as we already set some labels
    
    ax4.set_xlabel('Radial Distance (r/R)', fontsize=11, fontweight='bold')
    ax5.set_xlabel('Radial Distance (r/R)', fontsize=11, fontweight='bold')
    ax6.set_xlabel('Radial Distance (r/R)', fontsize=11, fontweight='bold')
    
    return fig

# Main execution
if __name__ == "__main__":
    print("Generating detailed annotated induction plots...\n")
    
    # Propeller configuration
    B = 3
    R = 0.5
    R_hub = 0.05
    rho = 1.225
    N = 50
    
    chord_fun = lambda r: 0.04
    theta_fun = lambda r: math.radians(20.0) * (1.0 - (r - R_hub)/(R - R_hub))
    
    # Generate plots for three conditions
    conditions = [
        (5.0, 4000, 'Low speed'),
        (10.0, 4000, 'Medium speed'),
        (15.0, 4000, 'High speed')
    ]
    
    for V, rpm, label in conditions:
        print(f"Generating: {label} (V={V} m/s, RPM={rpm})")
        fig = create_annotated_single_plot(B, chord_fun, theta_fun, R, R_hub, N, rho, V, rpm)
        
        filename = f"induction_detailed_V{V:.0f}_RPM{rpm}.png"
        fig.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"  ✓ Saved: {filename}\n")
    
    print("All detailed plots generated successfully!")
    plt.show()
