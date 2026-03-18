"""
This contains the implementation of scalar fields.
"""

import sympy as sp
from .base import Field
from ..spacetimes.base import Spacetime
from ..utils.tensor_op import covariant_derivative

class ScalarField(Field):
    '''
    Scalar field is a field that assigns a single scalar value to every point in spacetime.
    '''
    def __init__(self, spacetime, name="phi", mass=0, charge=0, xi=0):
        '''
        Parameters:
        spacetime: An instance of the Spacetime class that defines the geometry of the spacetime.
        name: The name of the scalar field (default is "phi").
        mass: The mass of the scalar field (default is 0, which means massless).
        charge: The charge of the scalar field (default is 0, which means uncharged).
        xi: The coupling constant to curvature (default is 0, which means minimal coupling).
        '''

        super().__init__(spacetime, name)
        self.field = sp.Function(name)(*spacetime.coords)  # φ(x^μ)
        self.mass = mass
        self.charge = charge
        self.xi = xi

    def equation_of_motion(self):
        '''
        The equation of motion for a scalar field in curved spacetime is given by the Klein-Gordon equation:
           (□ - m^2 - ξ R) φ = 0
        where □ is the d'Alembertian operator, m is the mass of the field, ξ is the coupling constant to curvature, and R is the Ricci scalar.
        '''
        spacetime = self.spacetime
        g = spacetime.metric()
        g_inv = spacetime.inverse_metric()
        coords = spacetime.coords

        # Calculate the d'Alembertian operator □φ = g^{μν} ∇_μ ∇_ν φ
        cov_derivative_1 = covariant_derivative(self.field, spacetime, index_type='lower')
        cov_derivative_2 = covariant_derivative(cov_derivative_1, spacetime, index_type='lower')

        dAlembertian = sum(g_inv[mu, nu] * cov_derivative_2[nu][mu] for mu in range(len(coords)) for nu in range(len(coords)))

        # Calculate the Ricci scalar R
        R = spacetime.ricci_scalar()

        # Construct the equation of motion
        eom = sp.simplify(dAlembertian - self.mass**2 * self.field - self.xi * R * self.field)

        return eom
    
    def lagrangian_density(self):
        '''
        The Lagrangian density for a scalar field in curved spacetime is given by:
           L = -1/2 g^{μν} ∂_μ φ ∂_ν φ - 1/2 m^2 φ^2 - 1/2 ξ R φ^2
        where g^{μν} is the inverse metric, m is the mass of the field, ξ is the coupling constant to curvature, and R is the Ricci scalar.
        '''
        spacetime = self.spacetime
        g_inv = spacetime.inverse_metric()
        coords = spacetime.coords

        # Calculate the kinetic term -1/2 g^{μν} ∂_μ φ ∂_ν φ
        kinetic_term = -0.5 * sum(g_inv[mu, nu] * sp.diff(self.field, coords[mu]) * sp.diff(self.field, coords[nu]) for mu in range(len(coords)) for nu in range(len(coords)))

        # Calculate the mass term -1/2 m^2 φ^2
        mass_term = -0.5 * self.mass**2 * self.field**2

        # Calculate the curvature coupling term -1/2 ξ R φ^2
        R = spacetime.ricci_scalar()
        curvature_coupling_term = -0.5 * self.xi * R * self.field**2

        # Construct the Lagrangian density
        lagrangian = sp.simplify(kinetic_term + mass_term + curvature_coupling_term)

        return lagrangian
    
    def energy_momentum_tensor(self):
        '''
        The energy-momentum tensor for a scalar field in curved spacetime is given by:
           T_{μν} = ∂_μ φ ∂_ν φ - g_{μν} (1/2 g^{αβ} ∂_α φ ∂_β φ + 1/2 m^2 φ^2 + 1/2 ξ R φ^2) + ξ (G_{μν} φ^2 - ∇_μ ∇_ν φ^2 + g_{μν} □ φ^2)
        where g_{μν} is the metric, m is the mass of the field, ξ is the coupling constant to curvature, R is the Ricci scalar, and G_{μν} is the Einstein tensor.
        '''
        spacetime = self.spacetime
        g = spacetime.metric()
        g_inv = spacetime.inverse_metric()
        coords = spacetime.coords

        # Calculate the kinetic term ∂_μ φ ∂_ν φ
        kinetic_term = [[sp.diff(self.field, coords[mu]) * sp.diff(self.field, coords[nu]) for nu in range(len(coords))] for mu in range(len(coords))]

        # Calculate the potential term -g_{μν} (1/2 g^{αβ} ∂_α φ ∂_β φ + 1/2 m^2 φ^2 + 1/2 ξ R φ^2)
        kinetic_contraction = 0.5 * sum(g_inv[alpha, beta] * sp.diff(self.field, coords[alpha]) * sp.diff(self.field, coords[beta]) for alpha in range(len(coords)) for beta in range(len(coords)))
        mass_term = 0.5 * self.mass**2 * self.field**2
        curvature_coupling_term = 0.5 * self.xi * spacetime.ricci_scalar() * self.field**2
        potential_term = [[-g[mu, nu] * (kinetic_contraction + mass_term + curvature_coupling_term) for nu in range(len(coords))] for mu in range(len(coords))]

        # Calculate the curvature coupling term ξ (G_{μν} φ^2 - ∇_μ ∇_ν φ^2 + g_{μν} □ φ^2)
        G_mu_nu = spacetime.einstein_tensor()
        cov_derivative_phi_squared = covariant_derivative(self.field**2, spacetime, index_type='lower')
        curvature_coupling_term = [[self.xi * (G_mu_nu[mu, nu] * self.field**2 - cov_derivative_phi_squared[nu][mu] + g[mu, nu] * sum(g_inv[alpha, beta] * cov_derivative_phi_squared[beta][alpha] for alpha in range(len(coords)) for beta in range(len(coords)))) for nu in range(len(coords))] for mu in range(len(coords))]

        # Construct the energy-momentum tensor
        T_mu_nu = [[sp.simplify(kinetic_term[mu][nu] + potential_term[mu][nu] + curvature_coupling_term[mu][nu]) for nu in range(len(coords))] for mu in range(len(coords))]

        return T_mu_nu
    
    def hamiltonian_density(self):
        '''
        The Hamiltonian density for a scalar field in curved spacetime can be derived from the Lagrangian density using a Legendre transformation. However, for simplicity, we can express it in terms of the canonical momentum π = ∂L/∂(∂_0 φ) and the field φ itself. The Hamiltonian density is given by:
           H = π ∂_0 φ - L
        where L is the Lagrangian density.
        '''
        spacetime = self.spacetime
        g_inv = spacetime.inverse_metric()
        coords = spacetime.coords

        # Calculate the canonical momentum π = ∂L/∂(∂_0 φ)
        pi = sp.diff(self.lagrangian_density(), sp.diff(self.field, coords[0]))

        # Calculate the Hamiltonian density H = π ∂_0 φ - L
        hamiltonian = sp.simplify(pi * sp.diff(self.field, coords[0]) - self.lagrangian_density())

        return hamiltonian
    
    def stress_energy_tensor(self):
        '''
        Alias for energy_momentum_tensor for compatibility with different terminologies. I may use different
        terms in different contexts, but they refer to the same physical quantity.
        '''
        return self.energy_momentum_tensor()