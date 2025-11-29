import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

class PolarDB:
    """
    Database for storing and interpolating 360-degree polars for multiple airfoils & Reynolds numbers.
    Structure:
        self.data[airfoil][Re] = dataframe with columns ['alpha','CL','CD','CM']
    """
    def __init__(self):
        self.data = {}

    def add_polar(self, airfoil: str, Re: float, df: pd.DataFrame):
        """
        Add one 360-degree polar for a specific airfoil and Reynolds number.
        df columns: ['alpha','CL','CD','CM']
        """
        airfoil = airfoil.lower()

        if airfoil not in self.data:
            self.data[airfoil] = {}

        self.data[airfoil][Re] = df.sort_values("alpha")

    def get_available_Re(self, airfoil: str):
        return sorted(self.data[airfoil.lower()].keys())

    def lookup_coeffs(self, airfoil: str, alpha: float, Re: float):
        """
        Return CL, CD, CM at given alpha and Re using:
        1) interpolation in alpha from 360-degree polars
        2) interpolation in Reynolds number (log-scale)
        """
        airfoil = airfoil.lower()

        if airfoil not in self.data:
            raise ValueError(f"Airfoil '{airfoil}' not found in database")

        # get all available Re for this airfoil
        available_Re = self.get_available_Re(airfoil)

        # if exact Re exists → use directly
        if Re in available_Re:
            df = self.data[airfoil][Re]
            return self._interp_alpha(df, alpha)

        # otherwise find bracketing Re values
        Re_array = np.array(available_Re)

        if Re < Re_array[0]:
            Re_low = Re_high = Re_array[0]
        elif Re > Re_array[-1]:
            Re_low = Re_high = Re_array[-1]
        else:
            # find nearest two values
            idx = np.searchsorted(Re_array, Re)
            Re_low = Re_array[idx - 1]
            Re_high = Re_array[idx]

        # get polars at the two Re values
        df_low = self.data[airfoil][Re_low]
        df_high = self.data[airfoil][Re_high]

        # interpolate α → get CL/CD at both Re levels
        CL_low, CD_low, CM_low = self._interp_alpha(df_low, alpha)
        CL_high, CD_high, CM_high = self._interp_alpha(df_high, alpha)

        # log-scale interpolation in Reynolds number
        if Re_low == Re_high:
            return CL_low, CD_low, CM_low

        t = (np.log(Re) - np.log(Re_low)) / (np.log(Re_high) - np.log(Re_low))

        CL = (1 - t) * CL_low + t * CL_high
        CD = (1 - t) * CD_low + t * CD_high
        CM = (1 - t) * CM_low + t * CM_high

        return CL, CD, CM

    def _interp_alpha(self, df: pd.DataFrame, alpha: float):
        """
        Interpolate CL, CD, CM at any alpha.
        df should contain the full 0–360° polar.
        """
        alpha_mod = alpha % 360.0

        interp_CL = interp1d(df["alpha"], df["CL"], kind="linear", fill_value="extrapolate")
        interp_CD = interp1d(df["alpha"], df["CD"], kind="linear", fill_value="extrapolate")
        interp_CM = interp1d(df["alpha"], df["CM"], kind="linear", fill_value="extrapolate")

        return float(interp_CL(alpha_mod)), float(interp_CD(alpha_mod)), float(interp_CM(alpha_mod))
