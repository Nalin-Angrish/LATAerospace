import numpy as np
import math
from extraction import getvalue

def getValue_2d():
    cl,aoa=getvalue()
    m0,alpha_l0=np.polyfit(cl,aoa,1)
    return (1/m0),alpha_l0

print(getValue_2d())