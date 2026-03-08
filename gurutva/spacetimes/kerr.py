"""
This module defines the Kerr spacetime class.
The Kerr spacetime describes the gravitational field outside a rotating massive object such as a rotating black hole.
The metric tensor for Kerr spacetime in Boyer-Lindquist coordinates (t, r, θ, φ) is given by:
   g_{μν} = [[-(1 - 2GM*r/ρ^2), 0, 0, -2GMa*r*sin^2(θ)/ρ^2],
              [0, ρ^2/Δ, 0, 0],
              [0, 0, ρ^2, 0],
              [-2GMa*r*sin^2(θ)/ρ^2, 0, 0, (r^2 + a^2 + 2GMa^2*r*sin^2(θ)/ρ^2)*sin^2(θ)]]
where ρ^2 = r^2 + a^2 * cos^2(θ) and Δ = r^2 - 2GM*r + a^2, and a is the spin parameter of the rotating object.
"""

from .base import Spacetime
import sympy as sp
from ..constants import G, M_SUN

class Kerr(Spacetime):
    def __init__(self, mass=M_SUN, spin=0): # Default mass is set to the solar mass, and spin is set to zero for a non-rotating object
        # Define the coordinates for Kerr spacetime (t, r, θ, φ)
        coords = sp.symbols('t r theta phi')
        super().__init__(coords) # Initialize the base Spacetime class with the defined coordinates
        self.mass = mass
        self.spin = spin

    def metric(self):
        t, r, theta, phi = self.coords
        M = self.mass # Mass of the central object (e.g., black hole or planet)
        a = self.spin # Spin parameter of the rotating object
        rho_squared = r**2 + a**2 * sp.cos(theta)**2
        delta = r**2 - 2 * G * M * r + a**2
        
        g_tt = -(1 - 2 * G * M * r / rho_squared)
        g_rr = rho_squared / delta
        g_theta_theta = rho_squared
        g_phi_phi = (r**2 + a**2 + 2 * G * M * a**2 * r * sp.sin(theta)**2 / rho_squared) * sp.sin(theta)**2
        g_tphi = -2 * G * M * a * r * sp.sin(theta)**2 / rho_squared
        
        g = sp.Matrix([[g_tt, 0, 0, g_tphi],
                       [0, g_rr, 0, 0],
                       [0, 0, g_theta_theta, 0],
                       [g_tphi, 0, 0, g_phi_phi]])
        
        return g