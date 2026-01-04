# Project Structure

This repository contains tools for **Propulsion (BEMT/XFOIL-based propeller analysis)** and **Wing Aerodynamics (VLM/VSM/DEP studies)**.

## Contributors
Work was done on this project by 8 students of IIT Ropar:
- [Dharmansh Vyas](https://www.linkedin.com/in/dharmansh-vyas-babab8258/) (4th Year, Mechanical Engineering)
- [Eakamjit Singh](https://www.linkedin.com/in/eakamjit-singh/) (3rd Year, Mechanical Engineering)
- [Nalin Angrish](https://www.linkedin.com/in/nalin-angrish/) (3rd Year, Mechanical Engineering)
- [Ponnathavan S](https://www.linkedin.com/in/ponnathavans/) (3rd Year, Engineering Physics)
- [Amey Kemkar](https://www.linkedin.com/in/amey-kemkar-2ab222313/) (2nd Year, Mechanical Engineering)
- [Aryan Deshmukh](https://www.linkedin.com/in/aryan-deshmukh-ba3667321/) (2nd Year, Metallurgical and Materials Engineering)
- [Ishan Gangwani](https://www.linkedin.com/in/ishan-gangwani-519424323/) (2nd Year, Mechanical Engineering)
- [Kian Sparrow](https://www.linkedin.com/in/kian-sparrow-4a6361320/) (2nd Year, Engineering Physics)

---

## 📁 Root Directory


### Root Files
- **requirements.txt** - Python dependencies for the project  
- **README.md** - Project documentation  
- **Additional Files** - Organized files for the project

---

## 🚀 `propulsion/` - Propeller & BEMT Analysis

Tools for **Propeller Design, Analysis, and Optimization** using **Blade Element Momentum Theory (BEMT)** and **XFOIL**.


### Key Files
- **Airfoils/** - Airfoil coordinate data for propeller sections  
- **Propeller/** - Our Arbitrary propeller geometry and related data  
- **bemt_dynamic_airfoil_optimization.py** - BEMT-based Collective and RPM optimization for an arbitrary propeller (with data given in `Propeller/`) 
- **bemt_general_optimization.py** - BEMT-based Collective and RPM optimization for a propeller with a general airfoil (NACA4412)  
- **bemt_rmit.py** - Comparison of our algorithm to the performance observed in the Paper by RMIT
- **bemt_with_xfoil.py** - BEMT coupled with XFOIL for viscous airfoil data  
- **propeller_data.csv** - The propeller performance measurements as observed in the Paper by RMIT  
- **xfoil_final.py** - Python wrapper for running XFOIL  
- **xfoil.exe** - XFOIL executable (Windows)

---

## ✈️ `wing/` - Wing Aerodynamics & DEP Studies

Contains **2D/3D wing aerodynamic solvers**, **DEP (Distributed Electric Propulsion)** simulations, and **vortex methods**.


### Key Files
- **Airfoils/** - Airfoil coordinate files for wing analysis  
- **DEP_2D.py** - 2D DEP aerodynamic analysis  
- **DEP_cruise.py** - Cruise condition analysis for our chosen DEP configuration  
- **DEP_takeoff.py** - Takeoff condition analysis for our chosen DEP configuration  
- **DEP_paper_2D_flat.py** - In the paper by Spence, they mentioned a formula for a flat blown airfoil. This evaluates that formula to get the `Cl` according to his observations.
- **helper.py** - Shared utility function to read camber line from csv
- **VLM_2d.ipynb** - 2D Vortex Lattice Method notebook  
- **VSM_wing.py** - Full Vortex Sheet Method wing model  
- **VSM_wing_simplified.py** - A slightly simplified version of VSM solver  

---

## 🛩️ `wing_no_jet/` - Wing Analysis without Jet/Prop Effects

Wing aerodynamic analysis **Excluding jet/propeller slipstream effects** for baseline comparisons.


### Key Files
- **Airfoils/** - Airfoil datasets  
- **Gamma.py** - Circulation distribution calculation  
- **extraction.py** - Data extraction and post-processing  
- **helper.py** - Utility functions  

---

## 🛩️ `propeller_simulation.zip` - ANSYS simulation of our propeller

Simulation of the optimal takeoff conditions as described by the output of `propulsion/bemt_dynamic_airfoil_optimization.py`.

## 🛩️ `wing_simulation.zip` - ANSYS simulation of our wing without blown effects

Simulation of a 3D wing as described by the output of `wing_no_jet/Gamma.py`.

## 🛩️ `dep_simulation.zip` - ANSYS simulation of our wing with blown effects

Simulation of the optimal wing as described by the output of `wing/DEP_takeoff.py`.

---

## ⚙️ Setup & Usage (Basic)

```bash
pip install -r requirements.txt
cd propulsion
python bemt_general_optimization.py
```
This can be extrapolated to run all other files as well.
