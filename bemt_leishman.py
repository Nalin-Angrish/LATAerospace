import math 
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Creating Functions for inflow ratio , thrust coeff and total thrust

def cal_inflow(r , swirl_inflow, Cl_alpha, solidity, local_pitch):
    inflow_term1 = math.sqrt(math.pow(((solidity * Cl_alpha / 16) - (swirl_inflow/2)),2) - (solidity * Cl_alpha * local_pitch * r / 8))
    inflow_term2 = (solidity * Cl_alpha / 16) - (swirl_inflow/2)
    inflow = inflow_term1 - inflow_term2
    return inflow

def elemental_thrust_coeff(r, dr, local_pitch, Cl_alpha, solidity, inflow):
    term1 = solidity * Cl_alpha * dr /2
    term2 = (local_pitch*math.pow(r,2)) - (inflow*r)
    elemtental_thrust = term1 * term2
    return elemtental_thrust

def total_thrust_coeff(r, dr , elemental_thrust_coeff):
    N = 1/dr
    sum = 0
    for i in range(N):
        sum = sum + (elemental_thrust_coeff[i])
    return sum

if __name__ == '__main__':
    pass
    

