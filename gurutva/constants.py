"""
Physical constants used throughout gurutva.

All constants are expressed in SI units unless otherwise stated.
"""

import numpy as np

############################
# Fundamental constants
############################

# Speed of light (m/s)
c = 299792458.0

# Gravitational constant (m^3 kg^-1 s^-2)
G = 6.67430e-11

# Planck constant (J s)
h = 6.62607015e-34

# Reduced Planck constant (J s)
hbar = h / (2 * np.pi)

# Boltzmann constant (J/K)
k_B = 1.380649e-23

# Stefan–Boltzmann constant (W m^-2 K^-4)
sigma = 5.670374419e-8


############################
# Planck scale
############################

m_P = np.sqrt(hbar * c / G)       # Planck mass (kg)
l_P = np.sqrt(hbar * G / c**3)    # Planck length (m)
t_P = np.sqrt(hbar * G / c**5)    # Planck time (s)
E_P = m_P * c**2                  # Planck energy (J)


############################
# Astronomical constants
############################

AU = 1.495978707e11               # Astronomical unit (m)
PC = 3.085677581491367e16         # Parsec (m)
LY = 9.4607304725808e15           # Light year (m)

M_SUN = 1.98847e30                # Solar mass (kg)
R_SUN = 6.957e8                   # Solar radius (m)


############################
# Mathematical constants
############################

PI = np.pi
TWO_PI = 2 * np.pi