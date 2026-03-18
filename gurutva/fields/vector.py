'''
This module defines the VectorField class, which represents a vector field in a given spacetime.
The VectorField class inherits from the Field class and provides additional functionality specific 
to vector fields, such as computing the field strength tensor and its dual.
'''

import sympy as sp
from .base import Field
from ..spacetimes.base import Spacetime
from ..utils.tensor_op import covariant_derivative

class VectorField(Field):
    '''
    A vector field assigns a vector to every point in spacetime, and can be used to describe various physical 
quantities such as electromagnetic fields, fluid velocities, or gravitational fields. 
    '''
    def __init__(self, spacetime, name="A", mass=0, charge=0, xi=0):
        '''
        Parameters:
        spacetime: An instance of the Spacetime class that defines the geometry of the spacetime.
        name: The name of the vector field (default is "A").
        mass: The mass of the vector field (default is 0, which means massless).
        charge: The charge of the vector field (default is 0, which means uncharged).
        xi: The coupling constant to curvature (default is 0, which means minimal coupling).
        '''

        super().__init__(spacetime, name)
        self.field = sp.Matrix(sp.Function(name)(*spacetime.coords) for _ in range(len(spacetime.coords)))  # A_μ(x^ν)
        self.mass = mass
        self.charge = charge
        self.xi = xi

    # Additional methods specific to vector fields can be added here, such as computing the field strength tensor F_{μν} = ∂_μ A_ν - ∂_ν A_μ and its dual.

    def field_strength_tensor(self):
        '''
        Compute the field strength tensor F_{μν} = ∂_μ A_ν - ∂_ν A_μ for the vector field.
        This is a fundamental quantity in electromagnetism and other gauge theories, and it encodes the dynamics of the vector field.
        '''
        spacetime = self.spacetime
        coords = spacetime.coords

        F = sp.Matrix(len(coords), len(coords), lambda i, j: sp.diff(self.field[j], coords[i]) - sp.diff(self.field[i], coords[j]))

        return F
    
    def dual_field_strength_tensor(self):
        '''
        Compute the dual of the field strength tensor, which is defined as *F^{μν} = 1/2 ε^{μνρσ} F_{ρσ}, where ε^{μνρσ} is the Levi-Civita symbol.
        The dual field strength tensor is important in the study of electromagnetic duality and topological properties of gauge fields.
        '''
        spacetime = self.spacetime
        coords = spacetime.coords
        F = self.field_strength_tensor()

        # Compute the dual using the Levi-Civita symbol
        epsilon = sp.tensor.array.MutableDenseNDimArray(sp.zeros(len(coords), len(coords), len(coords), len(coords)))
        for i in range(len(coords)):
            for j in range(len(coords)):
                for k in range(len(coords)):
                    for l in range(len(coords)):
                        epsilon[i, j, k, l] = sp.sign(sp.Permutation(i, j, k, l))

        dual_F = sp.Matrix(len(coords), len(coords), lambda mu, nu: 0.5 * sum(epsilon[mu, nu, rho, sigma] * F[rho, sigma] for rho in range(len(coords)) for sigma in range(len(coords))))

        return dual_F
    
    def lagrangian_density(self):
        '''
        The Lagrangian density for a vector field in curved spacetime can be given by:
           L = -1/4 F_{μν} F^{μν} - 1/2 m^2 A_μ A^μ - 1/2 ξ R A_μ A^μ
        where F_{μν} is the field strength tensor, m is the mass of the field, ξ is the coupling constant to curvature, and R is the Ricci scalar.
        '''
        spacetime = self.spacetime
        g_inv = spacetime.inverse_metric()
        coords = spacetime.coords

        F = self.field_strength_tensor()

        # Calculate the kinetic term -1/4 F_{μν} F^{μν}
        kinetic_term = -0.25 * sum(F[mu, nu] * sum(g_inv[mu, rho] * g_inv[nu, sigma] * F[rho, sigma] for rho in range(len(coords)) for sigma in range(len(coords))) for mu in range(len(coords)) for nu in range(len(coords)))

        # Calculate the mass term -1/2 m^2 A_μ A^μ
        mass_term = -0.5 * self.mass**2 * sum(g_inv[mu, nu] * self.field[mu] * self.field[nu] for mu in range(len(coords)) for nu in range(len(coords)))

        # Calculate the curvature coupling term -1/2 ξ R A_μ A^μ
        R = spacetime.ricci_scalar()
        curvature_coupling_term = -0.5 * self.xi * R * sum(g_inv[mu, nu] * self.field[mu] * self.field[nu] for mu in range(len(coords)) for nu in range(len(coords)))

        lagrangian = kinetic_term + mass_term + curvature_coupling_term

        return lagrangian
    
    def equation_of_motion(self):
        '''
        The equation of motion for a vector field in curved spacetime can be derived from the Lagrangian density using the Euler-Lagrange equations. 
        It can be expressed as:
           ∇_μ F^{μν} + m^2 A^ν + ξ R A^ν = 0
        where F^{μν} is the field strength tensor, m is the mass of the field, ξ is the coupling constant to curvature, and R is the Ricci scalar.
        '''
        spacetime = self.spacetime
        g_inv = spacetime.inverse_metric()
        coords = spacetime.coords

        F = self.field_strength_tensor()

        # Calculate ∇_μ F^{μν}
        cov_derivative_F = [[covariant_derivative(F[mu, nu], spacetime, index_type='upper') for nu in range(len(coords))] for mu in range(len(coords))]
        divergence_F = [sum(cov_derivative_F[mu][nu] for mu in range(len(coords))) for nu in range(len(coords))]

        # Calculate m^2 A^ν
        mass_term = [self.mass**2 * sum(g_inv[nu, mu] * self.field[mu] for mu in range(len(coords))) for nu in range(len(coords))]

        # Calculate ξ R A^ν
        R = spacetime.ricci_scalar()
        curvature_coupling_term = [self.xi * R * sum(g_inv[nu, mu] * self.field[mu] for mu in range(len(coords))) for nu in range(len(coords))]

        # Construct the equation of motion
        eom = [sp.simplify(divergence_F[nu] + mass_term[nu] + curvature_coupling_term[nu]) for nu in range(len(coords))]

        return eom
    
    def energy_momentum_tensor(self):
        self.coords = self.spacetime.coords
        g = self.spacetime.metric()
        g_inv = self.spacetime.inverse_metric()
        coords = self.spacetime.coords

        F = self.field_strength_tensor()

        # The energy-momentum tensor for a vector field can be derived from the Lagrangian density and is given by:
        # T_{μν} = F_{μα} F_{ν}^{ α} - 1/4 g_{μν} F_{αβ} F^{αβ} + m^2 (A_μ A_ν - 1/2 g_{μν} A_α A^α) + ξ (G_{μν} A_α A^α - ∇_μ ∇_ν A_α A^α + g_{μν} □ A_α A^α)
        T = sp.Matrix(len(coords), len(coords), lambda mu, nu: sum(F[mu, alpha] * sum(g_inv[nu, beta] * F[alpha, beta] for beta in range(len(coords))) for alpha in range(len(coords))) - 0.25 * g[mu, nu] * sum(F[alpha, beta] * sum(g_inv[alpha, rho] * g_inv[beta, sigma] * F[rho, sigma] for rho in range(len(coords)) for sigma in range(len(coords))) for alpha in range(len(coords)) for beta in range(len(coords))) + self.mass**2 * (self.field[mu] * self.field[nu] - 0.5 * g[mu, nu] * sum(g_inv[alpha, beta] * self.field[alpha] * self.field[beta] for alpha in range(len(coords)) for beta in range(len(coords)))) + self.xi * (sp.simplify(spacetime.einstein_tensor()[mu, nu] * sum(g_inv[alpha, beta] * self.field[alpha] * self.field[beta] for alpha in range(len(coords)) for beta in range(len(coords))) - covariant_derivative(covariant_derivative(sum(g_inv[alpha, beta] * self.field[alpha] * self.field[beta] for alpha in range(len(coords)) for beta in range(len(coords))), spacetime, index_type='lower')[mu][nu]) + g[mu, nu] * covariant_derivative(covariant_derivative(sum(g_inv[alpha, beta] * self.field[alpha] * self.field[beta] for alpha in range(len(coords)) for beta in range(len(coords))), spacetime, index_type='lower')[rho][rho]))))

        return T
    
    def stress_energy_tensor(self):
        '''
        Alias for energy_momentum_tensor for compatibility with different terminologies. I may use different
        terms in different contexts, but they refer to the same physical quantity.
        '''
        return self.energy_momentum_tensor()