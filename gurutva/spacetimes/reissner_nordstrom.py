"""
This module defines the Reissner-Nordström spacetime class.
The Reissner-Nordström spacetime describes the gravitational field outside a charged, non-rotating massive object such as a charged black hole.
The metric tensor for Reissner-Nordström spacetime in Schwarzschild coordinates (t, r, θ, φ) is given by:
   g_{μν} = diag(-(1 - 2GM/r + Q^2/r^2), (1 - 2GM/r + Q^2/r^2)^(-1), r^2, r^2 sin^2(θ))
where G is the gravitational constant, M is the mass of the central object, Q is its electric charge, and r is the radial coordinate.
"""

from .base import Spacetime
import sympy as sp
from ..constants import G, M_SUN

class ReissnerNordstrom(Spacetime):
    def __init__(self, mass=M_SUN, charge=0): # Default mass is set to the solar mass, and charge is set to zero for a neutral object
        # Define the coordinates for Reissner-Nordström spacetime (t, r, θ, φ)
        coords = sp.symbols('t r theta phi')
        super().__init__(coords) # Initialize the base Spacetime class with the defined coordinates
        self.mass = mass
        self.charge = charge

    def metric(self):
        t, r, theta, phi = self.coords
        M = self.mass # Mass of the central object (e.g., black hole or planet)
        Q = self.charge # Electric charge of the central object
        g_tt = -(1 - 2 * G * M / r + Q**2 / r**2)
        g_rr = (1 - 2 * G * M / r + Q**2 / r**2)**(-1)
        g_theta_theta = r**2
        g_phi_phi = r**2 * sp.sin(theta)**2
        
        g = sp.diag(g_tt, g_rr, g_theta_theta, g_phi_phi)
        return g