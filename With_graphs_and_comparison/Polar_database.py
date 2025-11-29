from xfoil_to_360 import XfoilPolar360
from polar_lookup import PolarDB
import pandas as pd
xfoil_path = r"C:\Users\ISHAN\Downloads\XFOIL6.99\xfoil.exe"
x360 = XfoilPolar360(xfoil_path=xfoil_path)

Re_list = [50000, 100000, 200000, 400000, 800000]
airfoil = "NACA2412"

pdb = PolarDB()

for Re in Re_list:
    out = x360.make_polar_and_extrapolate(
        airfoil,
        alphas=list(range(-10, 16)),
        Re=Re,
        Mach=0.0,
        AR=8,
        polar_filename=f"polar_Re{Re}.txt",
        out_csv=f"polar_Re{Re}_360.csv",
        plot_fig=None
    )

    df = pd.read_csv(f"polar_Re{Re}_360.csv")
    pdb.add_polar(airfoil, Re, df)

# save database for later use
import pickle
pickle.dump(pdb, open("polar_database.pkl", "wb"))
