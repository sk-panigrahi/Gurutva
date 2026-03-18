'''
This is the fields module of the library.
'''

class Field:
    def __init__(self, spacetime, name="F"):
        self.spacetime = spacetime
        self.coords = spacetime.coords
        self.name = name

        # Other properties
        self.mass = None
        self.charge = None
        self.spin = None
        self.xi = None  # Coupling constant to curvature

    def set_mass(self, mass):
        self.mass = mass
    
    def set_charge(self, charge):
        self.charge = charge

    def set_spin(self, spin):
        self.spin = spin

    def set_coupling_constant(self, xi):
        self.xi = xi

    # Additional stuff for field dynamics, interactions, etc.

    def equation_of_motion(self):
        '''
        Calculate the equation of motion for the field.
        This method should be implemented by subclasses to specify the dynamics of the field.
        '''
        raise NotImplementedError("Subclasses must implement the equation_of_motion method.")
    
    def energy_momentum_tensor(self):
        '''
        Calculate the energy-momentum tensor for the field.
        This method should be implemented by subclasses to specify the energy-momentum tensor of the field.
        '''
        raise NotImplementedError("Subclasses must implement the energy_momentum_tensor method.")
    
    def stress_energy_tensor(self):
        '''
        Alias for energy_momentum_tensor for compatibility with different terminologies. I may use different
        terms in different contexts, but they refer to the same physical quantity.
        '''
        return self.energy_momentum_tensor()
    
    def lagrangian_density(self):
        '''
        Calculate the Lagrangian density for the field.
        This method should be implemented by subclasses to specify the Lagrangian density of the field.
        '''
        raise NotImplementedError("Subclasses must implement the lagrangian_density method.")
    
    def hamiltonian_density(self):
        '''
        Calculate the Hamiltonian density for the field.
        This method should be implemented by subclasses to specify the Hamiltonian density of the field.
        '''
        raise NotImplementedError("Subclasses must implement the hamiltonian_density method.")
    