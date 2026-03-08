"""
This module defines the Schwarzschild spacetime class.
The Schwarzschild spacetime describes the gravitational field outside a spherical, 
non-rotating mass such as a black hole or a planet. The metric tensor for 
Schwarzschild spacetime in Schwarzschild coordinates (t, r, θ, φ) is given by:
   g_{μν} = diag(-(1 - 2GM/r), (1 - 2GM/r)^(-1), r^2, r^2 sin^2(θ))
"""

from .base import Spacetime
import sympy as sp
from ..constants import G, M_SUN

class Schwarzschild(Spacetime):
    def __init__(self, mass=M_SUN): # Default mass is set to the solar mass, but can be specified for other objects
        # Define the coordinates for Schwarzschild spacetime (t, r, θ, φ)
        coords = sp.symbols('t r theta phi')
        super().__init__(coords) # Initialize the base Spacetime class with the defined coordinates
        self.mass = mass

    def metric(self):
        t, r, theta, phi = self.coords
        M = self.mass # Mass of the central object (e.g., black hole or planet)
        g_tt = -(1 - 2 * G * M / r)
        g_rr = (1 - 2 * G * M / r)**(-1)
        g_theta_theta = r**2
        g_phi_phi = r**2 * sp.sin(theta)**2
        
        g = sp.diag(g_tt, g_rr, g_theta_theta, g_phi_phi)
        return g
    

