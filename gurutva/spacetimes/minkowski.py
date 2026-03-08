"""
This module defines the Minkowski spacetime class. 
The Minkowski spacetime is the simplest flat spacetime used in special relativity.
It is characterized by a metric tensor with signature (-+++), which describes an 
inertial frame of reference in the absence of gravity. The metric tensor for 
Minkowski spacetime is given by:
   g_{μν} = diag(-1, 1, 1, 1)
"""

from .base import Spacetime
import sympy as sp
from ..constants import G, M_SUN

class Minkowski(Spacetime):
    def __init__(self):
        # Define the coordinates for Minkowski spacetime (t, x, y, z)
        coords = sp.symbols('t x y z')
        super().__init__(coords) # Initialize the base Spacetime class with the defined coordinates

    def metric(self):
        # Define the Minkowski metric tensor g_{μν} = diag(-1, 1, 1, 1)
        g = sp.diag(-1, 1, 1, 1)
        return g