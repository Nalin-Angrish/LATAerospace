import math 
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# 1st make a function of transendental equation

def g(phi, solidity, V, w, r, cd_alpha, cl_alpha, local_pitch):
    h = math.sin(phi) + (0.25 * solidity * cd_alpha *(local_pitch - phi))
    e = 0.25 * solidity * cl_alpha * (local_pitch - phi)
    g = ((h*math.sin(phi)) - (e*math.cos(phi))) - (V / (w * r))*(h*math.cos(phi) + e*math.sin(phi))
    return g

# Next find its solution using regula falsi iteration
def iter_g(N, solidity, V, w, r, cd_alpha, cl_alpha, local_pitch):
    tol = 1e-6
    phi_left = 0
    phi_right = 0.5 * math.pi
    g_left = g(phi_left, solidity, V, w, r, cd_alpha, cl_alpha, local_pitch)
    g_right = g(phi_right, solidity, V, w, r, cd_alpha, cl_alpha, local_pitch)
    if g_left * g_right > 0:
        print("cant approach ahead")
        return
    else:
        for i in range(N):
            phi_new = phi_left - ((phi_right - phi_left)*g_left/(g_right - g_left))
            g_new = g(phi_new, solidity, V, w, r, cd_alpha, cl_alpha, local_pitch)
            if abs(g_new) < tol:
                # print(f"phi_new_on {i}th iter",phi_new)
                return phi_new
            else:
                if g_new * g_left < 0 :
                    phi_right = phi_new
                    g_right = g_new
                else:
                    phi_left = phi_new
                    g_left = g_new
        # print(f"Phi for {r} :",phi_new)
        return phi_new 

# function to calculate elemental thrust
def thrust_power_elemental(B, c, Cl_alpha, Cd_alpha, rho, V, w, r, local_pitch, no_of_iterations, R):
    solidity = B * c/(math.pi * 2 * r)
    phi = iter_g(no_of_iterations, solidity, V, w, r, Cd_alpha, Cl_alpha, local_pitch)
    alpha = local_pitch - phi
    f = (math.sin(phi) * (math.sin(phi) + 0.25 * solidity * Cd_alpha*(alpha))) - 0.25*solidity*Cl_alpha*(alpha)*math.cos(phi)
    # g = math.cos(phi) * (math.sin(phi) + 0.25 * solidity * Cd_alpha(local_pitch - phi)) + 0.25*solidity*Cl_alpha(local_pitch - phi)*math.sin(phi)

    Vres = V * math.sin(phi)/f
    dT = 0.5 * B * c * rho * (Vres**2) * (Cl_alpha*alpha*math.cos(phi) - Cd_alpha*alpha*math.sin(phi))
    dQ = 0.5 * B * c * rho * (Vres**2) * (Cl_alpha*alpha*math.sin(phi) + Cd_alpha*alpha*math.cos(phi)) * r
    print(f"r={r:.3f} phi={phi:.3f} alpha={alpha:.3f} sigma={solidity:.4f} dT={dT:.3f}")

    return dT , dQ

# summation of elemental thrust to get total thrust
def total_thrust_power(B, c, Cl_alpha, Cd_alpha, rho, V, w, local_pitch, no_of_iterations, R, No_of_elements, R_hub):
    N = No_of_elements - 1
    dr = (R - R_hub)/ N
    total_thrust = 0
    total_power = 0
    for i in range(N):
        r_n = R_hub + i*dr
        dt , dq= thrust_power_elemental(B, c, Cl_alpha, Cd_alpha, rho, V, w, r_n, local_pitch, no_of_iterations, R)
        if i == 0 :
            total_power += dq/2 
            total_thrust += dt/2
        else:
            total_power += dq 
            total_thrust += dt
    return total_thrust , total_power

if __name__ == '__main__':
    rpm = 6000
    w = rpm * 2 * math.pi/60
    V = 80
    B = 3
    c = 0.2
    R = 1.0
    local_pitch = math.pi * 10 / 180
    n_iter = 30
    N = 40
    Cl_alpha = 5
    Cd_alpha = 0.01
    rho = 1.225
    R_hub = 0.2
    t , q = total_thrust_power(B, c, Cl_alpha, Cd_alpha, rho, V, w , local_pitch, n_iter, R, N, R_hub)
    # t , q = thrust_power_elemental(B, c, Cl_alpha, Cd_alpha, rho, V, w, 0.3, local_pitch, n_iter, R)
    
    print("thrust: ",t)
    print("power: ",q)

    




    

    





