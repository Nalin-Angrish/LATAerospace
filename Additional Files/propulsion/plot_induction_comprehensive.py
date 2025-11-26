"""
Induction Factors and Inflow Angle Visualization for BEMT

This script generates comprehensive plots showing:
1. Axial induction factor (a) vs radial distance
2. Tangential induction factor (a') vs radial distance  
3. Inflow angle (φ) vs radial distance
4. Comparison across different operating conditions
"""

import numpy as np
import matplotlib.pyplot as plt
from bemt_rmit_improved import integrate_blade
import math

def plot_single_condition(B, chord_fun, theta_fun, R, R_hub, N, rho, V, rpm, 
                         pitch_beta=0.0, alpha0_fun=None):
    """
    Create a single comprehensive plot for one operating condition.
    """
    # Run BEMT
    T, Q, diags, nodes = integrate_blade(
        B, chord_fun, theta_fun, R, R_hub, N, rho, V, rpm,
        pitch_beta=pitch_beta, alpha0_fun=alpha0_fun
    )
    
    # Extract data
    r_R = np.array([d['r']/R for d in diags if d])
    a = np.array([d['a'] for d in diags if d])
    b = np.array([d['b'] for d in diags if d])
    phi_rad = np.array([d['phi'] for d in diags if d])
    phi_deg = np.degrees(phi_rad)
    
    # Create figure
    fig = plt.figure(figsize=(16, 10))
    
    # Plot 1: Axial Induction Factor
    ax1 = plt.subplot(2, 3, 1)
    ax1.plot(r_R, a, 'b-o', linewidth=2.5, markersize=6)
    ax1.fill_between(r_R, a, alpha=0.2, color='blue')
    ax1.set_xlabel('Radial Distance r/R', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Axial Induction Factor (a)', fontsize=11, fontweight='bold')
    ax1.set_title(f'Axial Induction Factor\nV={V} m/s, RPM={rpm}', 
                  fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim([0, 1])
    
    # Plot 2: Tangential Induction Factor
    ax2 = plt.subplot(2, 3, 2)
    ax2.plot(r_R, b, 'r-s', linewidth=2.5, markersize=6)
    ax2.fill_between(r_R, b, alpha=0.2, color='red')
    ax2.set_xlabel('Radial Distance r/R', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Tangential Induction Factor (a\')', fontsize=11, fontweight='bold')
    ax2.set_title(f'Tangential Induction Factor\nV={V} m/s, RPM={rpm}', 
                  fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim([0, 1])
    
    # Plot 3: Inflow Angle
    ax3 = plt.subplot(2, 3, 3)
    ax3.plot(r_R, phi_deg, 'g-^', linewidth=2.5, markersize=6)
    ax3.fill_between(r_R, phi_deg, alpha=0.2, color='green')
    ax3.set_xlabel('Radial Distance r/R', fontsize=11, fontweight='bold')
    ax3.set_ylabel('Inflow Angle φ (degrees)', fontsize=11, fontweight='bold')
    ax3.set_title(f'Inflow Angle\nV={V} m/s, RPM={rpm}', 
                  fontsize=12, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    ax3.set_xlim([0, 1])
    
    # Plot 4: Combined induction factors
    ax4 = plt.subplot(2, 3, 4)
    ax4.plot(r_R, a, 'b-o', linewidth=2.5, markersize=5, label='Axial (a)')
    ax4.plot(r_R, b, 'r-s', linewidth=2.5, markersize=5, label='Tangential (a\')')
    ax4.set_xlabel('Radial Distance r/R', fontsize=11, fontweight='bold')
    ax4.set_ylabel('Induction Factors', fontsize=11, fontweight='bold')
    ax4.set_title('Induction Factors Comparison', fontsize=12, fontweight='bold')
    ax4.legend(fontsize=10, loc='best')
    ax4.grid(True, alpha=0.3)
    ax4.set_xlim([0, 1])
    
    # Plot 5: Twist angle (theta) vs r/R
    ax5 = plt.subplot(2, 3, 5)
    theta_values = np.degrees(np.array([theta_fun(r) for r in r_R * R]))
    ax5.plot(r_R, theta_values, 'm-d', linewidth=2.5, markersize=5)
    ax5.fill_between(r_R, theta_values, alpha=0.2, color='magenta')
    ax5.set_xlabel('Radial Distance r/R', fontsize=11, fontweight='bold')
    ax5.set_ylabel('Blade Pitch Angle θ (degrees)', fontsize=11, fontweight='bold')
    ax5.set_title('Blade Twist Distribution', fontsize=12, fontweight='bold')
    ax5.grid(True, alpha=0.3)
    ax5.set_xlim([0, 1])
    
    # Plot 6: Summary table
    ax6 = plt.subplot(2, 3, 6)
    ax6.axis('off')
    
    # Calculate summary statistics
    Omega = 2 * math.pi * rpm / 60
    V_tip = Omega * R
    J = V / (rpm/60 * 2*R)  # Advance ratio
    
    summary_text = f"""
OPERATING CONDITION SUMMARY
{'='*45}

Input Parameters:
  Velocity:           {V:.2f} m/s
  RPM:                {rpm:.0f} rpm
  Advance Ratio (J):  {J:.4f}
  Collective Pitch:   {np.degrees(pitch_beta):.2f}°

Blade Geometry:
  Number of Blades:   {B}
  Radius:             {R:.3f} m
  Hub Radius:         {R_hub:.3f} m
  
Performance Metrics:
  Thrust (T):         {T:.4f} N
  Torque (Q):         {Q:.4f} N·m
  Input Power:        {Q*Omega:.2f} W
  Output Power:       {T*V:.2f} W
  
Induction Factors (Root/Mid/Tip):
  a:     {a[0]:.4f} / {a[len(a)//2]:.4f} / {a[-1]:.4f}
  a':    {b[0]:.4f} / {b[len(b)//2]:.4f} / {b[-1]:.4f}
  
Inflow Angle (Root/Mid/Tip):
  φ (deg): {phi_deg[0]:.2f}° / {phi_deg[len(phi_deg)//2]:.2f}° / {phi_deg[-1]:.2f}°
"""
    
    ax6.text(0.05, 0.95, summary_text, transform=ax6.transAxes,
            fontsize=9, verticalalignment='top', family='monospace',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    
    plt.tight_layout()
    return fig

def plot_multi_condition_comparison(B, chord_fun, theta_fun, R, R_hub, N, rho, 
                                   V_list, rpm):
    """
    Create comparison plots for multiple velocity conditions.
    """
    # Run BEMT for each condition
    results = {}
    for V in V_list:
        T, Q, diags, nodes = integrate_blade(
            B, chord_fun, theta_fun, R, R_hub, N, rho, V, rpm
        )
        
        r_R = np.array([d['r']/R for d in diags if d])
        a = np.array([d['a'] for d in diags if d])
        b = np.array([d['b'] for d in diags if d])
        phi = np.array([d['phi'] for d in diags if d])
        
        results[V] = {
            'r_R': r_R,
            'a': a,
            'b': b,
            'phi': np.degrees(phi),
            'T': T,
            'Q': Q
        }
    
    # Create figure
    fig = plt.figure(figsize=(16, 10))
    
    colors = plt.cm.tab10(np.linspace(0, 1, len(V_list)))
    
    # Plot 1: Axial induction factor comparison
    ax1 = plt.subplot(2, 2, 1)
    for i, V in enumerate(sorted(V_list)):
        ax1.plot(results[V]['r_R'], results[V]['a'], 
                color=colors[i], linewidth=2.5, marker='o', 
                markersize=4, label=f'{V:.1f} m/s')
    ax1.set_xlabel('Radial Distance r/R', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Axial Induction Factor (a)', fontsize=11, fontweight='bold')
    ax1.set_title('Axial Induction Factor - Multi-Condition', fontsize=12, fontweight='bold')
    ax1.legend(fontsize=10, title='Velocity', loc='best')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim([0, 1])
    
    # Plot 2: Tangential induction factor comparison
    ax2 = plt.subplot(2, 2, 2)
    for i, V in enumerate(sorted(V_list)):
        ax2.plot(results[V]['r_R'], results[V]['b'], 
                color=colors[i], linewidth=2.5, marker='s', 
                markersize=4, label=f'{V:.1f} m/s')
    ax2.set_xlabel('Radial Distance r/R', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Tangential Induction Factor (a\')', fontsize=11, fontweight='bold')
    ax2.set_title('Tangential Induction Factor - Multi-Condition', fontsize=12, fontweight='bold')
    ax2.legend(fontsize=10, title='Velocity', loc='best')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim([0, 1])
    
    # Plot 3: Inflow angle comparison
    ax3 = plt.subplot(2, 2, 3)
    for i, V in enumerate(sorted(V_list)):
        ax3.plot(results[V]['r_R'], results[V]['phi'], 
                color=colors[i], linewidth=2.5, marker='^', 
                markersize=4, label=f'{V:.1f} m/s')
    ax3.set_xlabel('Radial Distance r/R', fontsize=11, fontweight='bold')
    ax3.set_ylabel('Inflow Angle φ (degrees)', fontsize=11, fontweight='bold')
    ax3.set_title('Inflow Angle - Multi-Condition', fontsize=12, fontweight='bold')
    ax3.legend(fontsize=10, title='Velocity', loc='best')
    ax3.grid(True, alpha=0.3)
    ax3.set_xlim([0, 1])
    
    # Plot 4: Summary table of key metrics
    ax4 = plt.subplot(2, 2, 4)
    ax4.axis('off')
    
    # Create table data
    table_data = []
    for V in sorted(V_list):
        Omega = 2 * math.pi * rpm / 60
        T = results[V]['T']
        Q = results[V]['Q']
        J = V / (rpm/60 * 2*R)
        P_in = Q * Omega
        P_out = T * V
        eta = 100 * P_out / P_in if P_in > 0 else 0
        
        table_data.append([
            f'{V:.1f}',
            f'{T:.4f}',
            f'{Q:.4f}',
            f'{J:.4f}'
        ])
    
    # Create table text
    table_text = """VELOCITY COMPARISON @ RPM=4000
═════════════════════════════════════════
V (m/s)  Thrust(N)  Torque(N·m)  J(J)
─────────────────────────────────────────
"""
    for row in table_data:
        table_text += f"{row[0]:>8}  {row[1]:>9}  {row[2]:>11}  {row[3]:>5}\n"
    
    ax4.text(0.1, 0.95, table_text, transform=ax4.transAxes,
            fontsize=9, verticalalignment='top', family='monospace',
            bbox=dict(boxstyle='round', facecolor='lightcyan', alpha=0.8))
    
    plt.tight_layout()
    return fig

if __name__ == "__main__":
    # Propeller configuration
    B = 3
    R = 0.5
    R_hub = 0.05
    rho = 1.225
    N = 50
    
    # Blade geometry
    chord_fun = lambda r: 0.04
    theta_fun = lambda r: math.radians(20.0) * (1.0 - (r - R_hub)/(R - R_hub))
    
    print("Generating plots...")
    print("\n" + "="*60)
    print("PLOT 1: Single Operating Condition (V=10 m/s, RPM=4000)")
    print("="*60)
    fig1 = plot_single_condition(B, chord_fun, theta_fun, R, R_hub, N, rho,
                                V=10.0, rpm=4000)
    fig1.savefig('induction_factors_single_condition.png', dpi=300, bbox_inches='tight')
    print("✓ Saved: induction_factors_single_condition.png")
    
    print("\n" + "="*60)
    print("PLOT 2: Multi-Condition Comparison")
    print("="*60)
    V_conditions = [3.0, 5.0, 7.5, 10.0, 12.5, 15.0]
    fig2 = plot_multi_condition_comparison(B, chord_fun, theta_fun, R, R_hub, N, rho,
                                          V_conditions, rpm=4000)
    fig2.savefig('induction_factors_multi_condition.png', dpi=300, bbox_inches='tight')
    print("✓ Saved: induction_factors_multi_condition.png")
    
    print("\n" + "="*60)
    print("PLOT 3: High-Speed Condition (V=15 m/s, RPM=5000)")
    print("="*60)
    fig3 = plot_single_condition(B, chord_fun, theta_fun, R, R_hub, N, rho,
                                V=15.0, rpm=5000)
    fig3.savefig('induction_factors_high_speed.png', dpi=300, bbox_inches='tight')
    print("✓ Saved: induction_factors_high_speed.png")
    
    print("\n" + "="*60)
    print("All plots generated successfully!")
    print("="*60)
    
    plt.show()
