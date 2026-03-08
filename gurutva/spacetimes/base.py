"""
This module defines the base class for spacetimes in the Gurutva library.
It provides methods to calculate the metric tensor, inverse metric, determinant, 
Christoffel symbols, Ricci tensor, and Ricci scalar for a given spacetime. 
Subclasses should implement the `metric` method to specify the metric tensor for 
their particular spacetime.
"""


import sympy as sp

class Spacetime:
    def __init__(self, coords):
        self.coords = coords

    def metric(self):
        raise NotImplementedError("Subclasses must implement the metric method.")
    
    def inverse_metric(self):
        """
        Calculate the inverse of the metric tensor.
        The inverse metric g^{ij} is defined such that g^{ik} * g_{kj} = δ^i_j, 
        where δ^i_j is the Kronecker delta.
        """
        g = self.metric()
        return sp.simplify(g.inv())
    
    def determinant(self):
        g = self.metric()
        return sp.simplify(g.det())
    
    def christoffel_symbols(self):
        """
        Calculate the Christoffel symbols of the second kind.
        The Christoffel symbols are given by:
           Γ^i_{jk} = (1/2) * g^{im} * (∂_j g_{mk} + ∂_k g_{mj} - ∂_ m g_{jk})
        where g^{im} is the inverse metric and g_{jk} is the metric tensor.
        """
        g = self.metric()
        g_inv = self.inverse_metric()
        coords = self.coords
        n = len(coords)
        
        Gamma = [[[0 for _ in range(n)] for _ in range(n)] for _ in range(n)]
        
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    Gamma[i][j][k] = sp.Rational(1, 2) * sum(
                        g_inv[i, m] * (sp.diff(g[m, j], coords[k]) + sp.diff(g[m, k], coords[j]) - sp.diff(g[j, k], coords[m]))
                        for m in range(n)
                    )
        
        return Gamma
    
    def ricci_tensor(self):
        """
        Calculate the Ricci tensor R_{ij} from the Christoffel symbols.
        The Ricci tensor is given by:
           R_{ij} = ∂_k Γ^k_{ij} - ∂_j Γ^k_{ik} + Γ^k_{ij} Γ^m_{km} - Γ^k_{ik} Γ^m_{jm}
        where Γ^k_{ij} are the Christoffel symbols.
        """
        Gamma = self.christoffel_symbols()
        coords = self.coords
        n = len(coords)
        
        R = [[0 for _ in range(n)] for _ in range(n)]
        
        for i in range(n):
            for j in range(n):
                R[i][j] = sum(
                    sp.diff(Gamma[k][i][j], coords[k]) - sp.diff(Gamma[k][i][k], coords[j]) +
                    sum(Gamma[k][i][m] * Gamma[m][j][k] - Gamma[k][i][k] * Gamma[m][j][m] for m in range(n))
                    for k in range(n)
                )
        
        return R
    
    def ricci_scalar(self):
        """
        Calculate the Ricci scalar R by contracting the Ricci tensor with the inverse metric.
        The Ricci scalar is given by:
           R = g^{ij} R_{ij}
        """
        R = self.ricci_tensor()
        g_inv = self.inverse_metric()
        return sp.simplify(sum(g_inv[i, j] * R[i][j] for i in range(len(self.coords)) for j in range(len(self.coords))))
