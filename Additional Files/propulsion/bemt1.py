"""
Simple script demonstrating an iterative solver for the axial (a) and
angular (a_dash) induction factors for a single blade element.

This implements a fixed-point iteration that updates a and a_dash while
recomputing the inflow angle phi, the Prandtl tip-loss factor F and the
force coefficients CT and CQ (built from Cl and Cd). The implementation
accepts either numeric Cl/Cd values or callables Cl(alpha)/Cd(alpha).

This file is intended as an educational / standalone demonstration and is
not integrated into the full pyBEMT library. It can be adapted to call
airfoil lookup functions where available.
"""

import math

from typing import Callable, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def prandtl_tip_loss(
    B: int, R: float, r: float, phi: float, radius_hub: float = 0.0
) -> float:
    """Prandtl tip+hub loss factor."""

    def prandtl_one(dr, rloc, phi_loc):
        if sin_phi == 0:
            return 1.0
        f = B * dr / (2.0 * rloc * sin_phi)
        if -f > 500:
            return 1.0
        return 2.0 * math.acos(min(1.0, math.exp(-f))) / math.pi

    sin_phi = math.sin(phi)
    if phi == 0:
        return 1.0
    Ftip = prandtl_one(R - r, r, phi)
    Fhub = prandtl_one(r - radius_hub, r, phi)
    return Ftip * Fhub


def solve_induction(
    V: float,
    omega: float,
    r: float,
    chord: float,
    B: int,
    pitch: float,
    Cl: "Callable[[float], float] | float",
    Cd: "Callable[[float], float] | float",
    R: float,
    radius_hub: float = 0.0,
    rho: float = 1.225,
    tol: float = 1e-6,
    maxiter: int = 200,
    relax: float = 0.25,
) -> Tuple[float, float, float]:
    """Iteratively solve for axial a, angular a_dash and inflow angle phi.

    Inputs:
      V, omega: inflow and rotational velocities
      r: local radius
      chord: blade chord at the section
      B: blade count
      pitch: blade pitch (rad)
      Cl, Cd: either callables Cl(alpha) and Cd(alpha) or float constants
      R: rotor radius (for tip loss)

    Returns (a, a_dash, phi)
    """
    # solidity
    sigma = B * chord / (2.0 * math.pi * r)

    # helper to evaluate Cl/Cd whether they are callables or constants
    def eval_Cl(alpha):
        return Cl(alpha) if callable(Cl) else Cl

    def eval_Cd(alpha):
        return Cd(alpha) if callable(Cd) else Cd

    # initial guesses
    a = 0.0
    ap = 0.0

    for _ in range(maxiter):
        # compute inflow angle from current a, ap
        denom = omega * r * (1 - ap)
        # if denom is zero (stalled rotation), avoid division by zero
        if denom == 0.0:
            phi = 0.5 * math.pi
        else:
            phi = math.atan((1 + a) * V / denom)

        # Tip/hub loss
        F = prandtl_tip_loss(B, R, r, phi, radius_hub)

        # aerodynamic coefficients at current angle of attack
        alpha = pitch - phi
        cl = eval_Cl(alpha)
        cd = eval_Cd(alpha)

        # force coefficients CT and CQ (axial and tangential components)
        CT = cl * math.cos(phi) - cd * math.sin(phi)
        CQ = cl * math.sin(phi) + cd * math.cos(phi)

        sin_phi = math.sin(phi)
        cos_phi = math.cos(phi)

        # avoid division by zero; when phi->0 use algebraic limits similar to library
        if abs(sin_phi) < 1e-8:
            # use algebraic fallback consistent with rotor.py
            # kappa = 4*F*sin^2/(sigma*CT) -> a = 1/(kappa - 1) with C=1 (rotor)
            if CT == 0.0:
                anew = a
            else:
                kappa = (
                    4.0 * F * sin_phi * sin_phi / (sigma * CT)
                    if sigma * CT != 0
                    else float("inf")
                )
                anew = 1.0 / (kappa - 1.0)
        else:
            RHS_a = (sigma * CT * F) / (4.0 * sin_phi * sin_phi)
            # fixed point closed form
            if RHS_a <= -1e-12:
                # guard against division by (1+RHS_a) ~ 0
                anew = a
            else:
                anew = RHS_a / (1.0 + RHS_a)

        # a' update
        if abs(sin_phi * cos_phi) < 1e-12:
            # fallback using library expression
            if CQ == 0.0:
                apnew = ap
            else:
                kappap = (
                    4.0 * F * sin_phi * cos_phi / (sigma * CQ)
                    if sigma * CQ != 0
                    else float("inf")
                )
                apnew = 1.0 / (kappap + 1.0)
        else:
            RHS_ap = (sigma * CQ * F) / (4.0 * sin_phi * cos_phi)
            # solve a'/(1 + a') = RHS_ap  => a' = RHS_ap/(1 - RHS_ap)
            if abs(1.0 - RHS_ap) < 1e-12:
                apnew = ap
            else:
                apnew = RHS_ap / (1.0 - RHS_ap)

        # relaxation
        a = relax * anew + (1.0 - relax) * a
        ap = relax * apnew + (1.0 - relax) * ap

        # convergence
        if abs(a - anew) < tol and abs(ap - apnew) < tol:
            break

    return a, ap, phi


if __name__ == "__main__":
    # Example usage with constant Cl/Cd
    V = 10.0  # m/s
    rpm = 200.0
    omega = rpm * 2.0 * math.pi / 60.0
    r = 0.15
    chord = 0.2
    B = 3
    pitch = 3.0
    R = 1.0

    pathname_lift = "../data/naca2412_lift_coeff.csv"
    pathname_drag = "../data/naca2412_drag_coeff.csv"
    # simple constant lift/drag for demo
    try:
        Cl_val = pd.read_csv(pathname_lift, header=None, skiprows=1)
        x = 1/0
        Cd_val = pd.read_csv(pathname_drag, header=None, skiprows=1)
    except Exception as e:
        print(f"Error reading airfoil data: {e}")
        Cl_val = 0.5  # fallback constant
        Cd_val = 0.01  # fallback constant

    def element_loading(
        V, omega, r, chord, B, pitch, Cl, Cd, R, radius_hub=0.0, rho=1.225
    ):
        """Compute dT/dr and dQ/dr at radius r (spanwise loading)."""
        a, ap, phi = solve_induction(
            V, omega, r, chord, B, pitch, Cl, Cd, R, radius_hub, rho
        )

        # aerodynamic coefficients at current AoA
        alpha = pitch - phi
        cl = Cl(alpha) if callable(Cl) else Cl
        cd = Cd(alpha) if callable(Cd) else Cd

        CT = cl * math.cos(phi) - cd * math.sin(phi)
        CQ = cl * math.sin(phi) + cd * math.cos(phi)

        v_axial = (1 + a) * V
        v_tang = (1 - ap) * omega * r
        U = math.hypot(v_axial, v_tang)

        # blade element theory (force from one blade times B blades)
        dTdr = 0.5 * rho * U**2 * chord * CT * B
        dQdr = 0.5 * rho * U**2 * chord * CQ * B * r

        return dTdr, CT, a, ap, phi

    # Generate spanwise loading plots
    r_hub = 0.02  # hub radius
    r = np.linspace(r_hub, R, 100)  # radial positions to evaluate
    dTdr = np.zeros_like(r)
    CT = np.zeros_like(r)
    a_dist = np.zeros_like(r)
    ap_dist = np.zeros_like(r)
    phi_dist = np.zeros_like(r)

    for i, ri in enumerate(r):
        dTdr[i], CT[i], a_dist[i], ap_dist[i], phi_dist[i] = element_loading(
            V, omega, ri, chord, B, pitch, Cl_val, Cd_val, R, r_hub
        )

    # Plot spanwise distributions
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle("Spanwise Distributions")

    # Loading distributions
    ax1.plot(r / R, dTdr, "b-", label="dT/dr")
    ax1.set_xlabel("r/R")
    ax1.set_ylabel("dT/dr (N/m)")
    ax1.grid(True)
    ax1.legend()

    ax2.plot(r / R, CT, "r-", label="dQ/dr")
    ax2.set_xlabel("r/R")
    ax2.set_ylabel("dQ/dr (Nm/m)")
    ax2.grid(True)
    ax2.legend()

    # Induction factors
    ax3.plot(r / R, a_dist, "g-", label="a")
    ax3.plot(r / R, ap_dist, "b--", label="a'")
    ax3.set_xlabel("r/R")
    ax3.set_ylabel("Induction Factors")
    ax3.grid(True)
    ax3.legend()

    # Flow angle
    ax4.plot(r / R, np.degrees(phi_dist), "k-", label="φ")
    ax4.set_xlabel("r/R")
    ax4.set_ylabel("Flow Angle φ (deg)")
    ax4.grid(True)
    ax4.legend()

    plt.tight_layout()
    plt.show()

    # Print some key values
    print("\nSpanwise loading at r/R = 0.75:")
    idx = np.argmin(np.abs(r / R - 0.75))
    print(f"dT/dr = {dTdr[idx]:.2f} N/m")
    # print(f'dQ/dr = {dQdr[idx]:.2f} Nm/m')
    print(f"a = {a_dist[idx]:.3f}")
    print(f"a' = {ap_dist[idx]:.3f}")
    print(f"phi = {np.degrees(phi_dist[idx]):.1f} deg")

    # Compute total thrust, torque and power by numerical integration
    dr = r[1] - r[0]  # spacing between points
    T = np.sum(dTdr * dr)  # integrate dT/dr
    # Q = np.sum(dQdr * dr)  # integrate dQ/dr
    # P = Q * omega

    print("\nTotal integrated values:")
    print(f"Thrust (N): {T:.2f}")
    # print(f'Torque (Nm): {Q:.2f}')
    # print(f'Power (W): {P:.2f}')

    # Save results to CSV
    results = pd.DataFrame(
        {
            "r/R": r / R,
            "r (m)": r,
            "dT/dr (N/m)": dTdr,
            # 'dQ/dr (Nm/m)': dQdr,
            "a": a_dist,
            "a_prime": ap_dist,
            "phi (deg)": np.degrees(phi_dist),
        }
    )

    csv_filename = "../out/blade_loading.csv"
    results.to_csv(csv_filename, index=False)
    print(f"\nSpanwise distributions saved to {csv_filename}")
