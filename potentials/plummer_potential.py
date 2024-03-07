import numpy as np
from amuse.lab import *

class PlummerPotential(object):
    def __init__(self):
        self.mass= 0 | units.MSun
        self.radius = 0 | units.pc
        self.epsilon2 = (0 | units.parsec)**2

    def set_parameters(self, mass, radius, eps2=(0|units.pc)**2):
        self.mass = mass
        self.radius = radius
        self.epsilon2 = eps2
        
    def get_potential_at_point(self,eps,x,y,z):
        r=(x**2+y**2+z**2 + self.radius**2 + self.epsilon2)**0.5
        potential = 2 * constants.G * self.mass/r
        return potential

    def get_gravity_at_point(self, eps, x,y,z):
        phi_0 = self.get_potential_at_point(eps, x,y,z)
        #grav = AdaptingVectorQuantity()
        dpos = 0.001*(x**2+y**2+z**2).sqrt()
        phi_dx = self.get_potential_at_point(0,x+dpos,y,z) - phi_0
        phi_dy = self.get_potential_at_point(0,x,y+dpos,z) - phi_0
        phi_dz = self.get_potential_at_point(0,x,y, z+dpos) - phi_0
        return phi_dx/dpos, phi_dy/dpos, phi_dz/dpos

    def mass_in(self, r):
        return self.mass * r**3/(r**2 + self.radius**2)**(3./2.)

    def density(self, r):
        return 3.*self.mass/(4*np.pi*self.radius**3) * (1 + (r/self.radius)**2)**(-5./2.)

    def velocity_dispersion(self, R):
        r=(R**2 + self.radius**2)**0.5
        return np.sqrt(constants.G*self.mass_in(r)/(6*r))

    def circular_velocity(self, R):
        return np.sqrt(constants.G*self.mass_in(R)/R)
    
    def stop(self):
        return
