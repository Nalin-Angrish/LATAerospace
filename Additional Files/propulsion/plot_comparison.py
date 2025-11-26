import numpy as np
import matplotlib.pyplot as plt
from bemt_rmit_improved import integrate_blade
import math

def plot_comparison_conditions(B, chord_fun, theta_fun, R, R_hub, N, rho, 
                               V_conditions, rpm, pitch_beta=0.0):
    """
    Plot induction factors and inflow angle for multiple velocity conditions.
    
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
    V_conditions : list
        List of freestream velocities to compare (m/s)
    rpm : float
        Rotational speed (rpm)
    pitch_beta : float
        Collective pitch angle (rad)
    """
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Colors for different conditions
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    
    # Run BEMT for each velocity condition
    all_data = {}
    for i, V in enumerate(V_conditions):
        T, Q, diags, nodes = integrate_blade(
            B, chord_fun, theta_fun, R, R_hub, N, rho, V, rpm,
            pitch_beta=pitch_beta
        )
        
        r_R = np.array([d['r']/R for d in diags if d])
        a = np.array([d['a'] for d in diags if d])
        b = np.array([d['b'] for d in diags if d])
        phi = np.array([d['phi'] for d in diags if d])
        
        all_data[V] = {
            'r_R': r_R,
            'a': a,
            'b': b,
            'phi': np.degrees(phi),
            'T': T,
            'Q': Q
        }
    
    # Plot 1: Axial Induction Factor
    ax1 = axes[0, 0]
    for i, V in enumerate(V_conditions):
        ax1.plot(all_data[V]['r_R'], all_data[V]['a'], 
                color=colors[i % len(colors)], linewidth=2.5, marker='o', 
                markersize=4, label=f'V={V} m/s')
    ax1.set_xlabel('Radial Distance (r/R)', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Axial Induction Factor (a)', fontsize=11, fontweight='bold')
    ax1.set_title('Axial Induction Factor vs Radial Distance', fontsize=12, fontweight='bold')
    ax1.legend(fontsize=10, loc='best')
    ax1.grid(True, alpha=0.3, linestyle='--')
    ax1.set_xlim([0, 1])
    
    # Plot 2: Tangential Induction Factor
    ax2 = axes[0, 1]
    for i, V in enumerate(V_conditions):
        ax2.plot(all_data[V]['r_R'], all_data[V]['b'], 
                color=colors[i % len(colors)], linewidth=2.5, marker='s', 
                markersize=4, label=f'V={V} m/s')
    ax2.set_xlabel('Radial Distance (r/R)', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Tangential Induction Factor (a\')', fontsize=11, fontweight='bold')
    ax2.set_title('Tangential Induction Factor vs Radial Distance', fontsize=12, fontweight='bold')
    ax2.legend(fontsize=10, loc='best')
    ax2.grid(True, alpha=0.3, linestyle='--')
    ax2.set_xlim([0, 1])
    
    # Plot 3: Inflow Angle
    ax3 = axes[1, 0]
    for i, V in enumerate(V_conditions):
        ax3.plot(all_data[V]['r_R'], all_data[V]['phi'], 
                color=colors[i % len(colors)], linewidth=2.5, marker='^', 
                markersize=4, label=f'V={V} m/s')
    ax3.set_xlabel('Radial Distance (r/R)', fontsize=11, fontweight='bold')
    ax3.set_ylabel('Inflow Angle (degrees)', fontsize=11, fontweight='bold')
    ax3.set_title('Inflow Angle (φ) vs Radial Distance', fontsize=12, fontweight='bold')
    ax3.legend(fontsize=10, loc='best')
    ax3.grid(True, alpha=0.3, linestyle='--')
    ax3.set_xlim([0, 1])
    
    # Plot 4: Performance comparison
    ax4 = axes[1, 1]
    V_list = sorted(list(all_data.keys()))
    T_list = [all_data[V]['T'] for V in V_list]
    P_in_list = [all_data[V]['Q'] * 2*math.pi*rpm/60 for V in V_list]
    P_out_list = [all_data[V]['T'] * V for V in V_list]
    
    ax4_thrust = ax4
    ax4_power = ax4.twinx()
    
    l1 = ax4_thrust.plot(V_list, T_list, 'b-o', linewidth=2.5, markersize=8, label='Thrust')
    ax4_thrust.set_xlabel('Freestream Velocity (m/s)', fontsize=11, fontweight='bold')
    ax4_thrust.set_ylabel('Thrust (N)', fontsize=11, fontweight='bold', color='blue')
    ax4_thrust.tick_params(axis='y', labelcolor='blue')
    
    l2 = ax4_power.plot(V_list, P_in_list, 'r--s', linewidth=2.5, markersize=8, label='Input Power')
    l3 = ax4_power.plot(V_list, P_out_list, 'g--^', linewidth=2.5, markersize=8, label='Output Power')
    ax4_power.set_ylabel('Power (W)', fontsize=11, fontweight='bold', color='red')
    ax4_power.tick_params(axis='y', labelcolor='red')
    
    ax4.set_title('Performance vs Velocity (RPM=4000)', fontsize=12, fontweight='bold')
    ax4.grid(True, alpha=0.3, linestyle='--')
    
    # Add combined legend
    lines = l1 + l2 + l3
    labels = [l.get_label() for l in lines]
    ax4.legend(lines, labels, fontsize=10, loc='upper left')
    
    plt.tight_layout()
    return fig, all_data

if __name__ == "__main__":
    # Propeller parameters
    B = 3
    R = 0.5
    R_hub = 0.05
    rho = 1.225
    N = 40
    
    chord_fun = lambda r: 0.04
    theta_fun = lambda r: math.radians(20.0) * (1 - (r - R_hub)/(R - R_hub))
    
    # Compare multiple velocity conditions at constant RPM
    print("Generating comparison plot for multiple velocity conditions...")
    V_conditions = [3.0, 5.0, 7.5, 10.0, 12.5, 15.0]  # m/s
    
    fig, data = plot_comparison_conditions(
        B, chord_fun, theta_fun, R, R_hub, N, rho,
        V_conditions, rpm=4000
    )
    
    fig.savefig('induction_factors_comparison.png', dpi=300, bbox_inches='tight')
    
    # Print summary table
    print("\n" + "="*80)
    print("INDUCTION FACTORS SUMMARY ACROSS VELOCITY CONDITIONS")
    print("="*80)
    print(f"{'Velocity':<12} {'Thrust':<12} {'Input Power':<15} {'Output Power':<15} {'Efficiency':<12}")
    print(f"{'(m/s)':<12} {'(N)':<12} {'(W)':<15} {'(W)':<15} {'(%)':<12}")
    print("-"*80)
    
    for V in sorted(V_conditions):
        T = data[V]['T']
        Q = data[V]['Q']
        Omega = 2*math.pi*4000/60
        P_in = Q * Omega
        P_out = T * V
        eta = 100 * P_out / P_in if P_in > 0 else 0
        
        print(f"{V:<12.1f} {T:<12.3f} {P_in:<15.2f} {P_out:<15.2f} {eta:<12.1f}")
    
    print("="*80)
    print("\nPlot saved as: induction_factors_comparison.png")
    plt.show()
