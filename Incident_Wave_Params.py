"""This file contains parameters for the incident wave used in scattering problems.
Note: We cannot just define these in Main.py for reasons to do with a feature of bempp-cl."""

import numpy as np

INC_WL = 1  # Wavelength of the incident plane wave in units of UNIT as defined in Main.py
INC_DIR = np.array([0, 0, 1.0])  # Direction of the incident plane wave. Must be a unit vector. Has no units.
INC_POL = np.array([1.0, 0, 0])  # Polarization of the incident plane wave. Must be a unit vector. Has no units.
