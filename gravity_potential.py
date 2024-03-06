import numpy as np
from amuse.datamodel import Particles
from amuse.units import units
from amuse.units import constants
from amuse.units import nbody_system
from amuse.couple import bridge
from amuse.community.ph4.interface import Ph4
from matplotlib import pyplot as plt
from amuse.ic.plummer import new_plummer_model
from amuse.ic.salpeter import new_salpeter_mass_distribution
from amuse.lab import read_set_from_file

from plot.plot_cluster import plot_cluster
from read_parameters import new_option_parser

class CompositeGravityCode(object):
    def __init__(self, converter):
        self.model_time = 0 | units.Myr
        self.particles = Particles()
        self.converter = converter
        self.potential = PointParticlePotential()
        #self.potential = PlummerPotential()

    def evolve_model(self, model_time, verbose=False):
        if verbose:
            print(f"evolve_plummer to t={model_time.in_(units.Myr)}: x={self.particles.x}")
        dt = model_time - self.model_time
        self.particles.position += self.particles.velocity*dt
        self.model_time = model_time

    def add_particle(self, particle):
        self.particles.add_particle(particle)
        self.potential.set_parameters(particle.mass, particle.radius)
        
    @property
    def potential_energy(self):
        return quantities.zero
    @property 
    def kinetic_energy(self):
        return (0.5*self.particles.mass \
                   *self.particles.velocity.lengths()**2).sum()

    def get_potential_at_point(self,eps,x,y,z):
        x = x - self.particles[0].x
        y = y - self.particles[0].y
        z = z - self.particles[0].z
        plummer_potential = self.potential.get_potential_at_point(eps,x,y,z)
        return plummer_potential

    def get_gravity_at_point(self, eps, x,y,z):
        x = x - self.particles[0].x
        y = y - self.particles[0].y
        z = z - self.particles[0].z
        phi_dx, phi_dy, phi_dz = self.potential.get_gravity_at_point(eps, x,y,z)
        return phi_dx, phi_dy, phi_dz

    def stop(self):
        self.potential.stop()                
        return

class PointParticlePotential(object):
    def __init__(self):
        self.mass= 0 | units.MSun
        self.radius = 0 | units.pc
        self.epsilon2 = (0 | units.parsec)**2

    def set_parameters(self, mass, radius):
        self.mass = mass
        self.radius = radius
        
    def get_potential_at_point(self,eps,x,y,z):
        r=(x**2+y**2+z**2 + self.epsilon2)**0.5
        potential = constants.G * self.mass/r
        return potential

    def get_gravity_at_point(self, eps, x,y,z):
        phi_0 = self.get_potential_at_point(eps, x,y,z)
        dpos = 0.001*(x**2+y**2+z**2).sqrt()
        phi_dx = self.get_potential_at_point(0,x+dpos,y,z) - phi_0
        phi_dy = self.get_potential_at_point(0,x,y+dpos,z) - phi_0
        phi_dz = self.get_potential_at_point(0,x,y, z+dpos) - phi_0
        return phi_dx/dpos, phi_dy/dpos, phi_dz/dpos

    def mass_in(self, r):
        return self.mass * r**3/(r**2 + self.radius**2)**(3./2.)

    def density(self, r):
        return 3.*self.mass/(4*np.pi*self.radius**2)

    def velocity_dispersion(self, R):
        return np.sqrt(constants.G*self.mass/R)

    def circular_velocity(self, R):
        return np.sqrt(constants.G*self.mass/R)
    
    def stop(self):
        return
    
class PlummerPotential(object):

    def __init__(self):
        self.mass= 0 | units.MSun
        self.radius = 0 | units.pc
        self.epsilon2 = (0 | units.parsec)**2

    def set_parameters(self, mass, radius):
        self.mass = mass
        self.radius = radius
        
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

if __name__ == "__main__":
    o, arguments = new_option_parser().parse_args()

    particles = read_set_from_file(o.filename, close_file=True)
    m = particles.mass.sum()
    r = particles.position.length().sum()/len(particles)
    converter=nbody_system.nbody_to_si(m, r)

    gravity = []
    channel = {"to_gr":[],
               "from_gr":[]}
    for pi in range(len(particles)):
        gravity.append(CompositeGravityCode(converter))    
        gravity[-1].add_particle(particles[pi].as_set())
        channel["to_gr"].append(particles.new_channel_to(gravity[-1].particles))
        channel["from_gr"].append(gravity[-1].particles.new_channel_to(particles))

        
    system=bridge.Bridge(verbose=False)
    for gi in range(len(gravity)):
        for gj in range(len(gravity)):
            if gi != gj:
                #print(f"Add code: ({gi}, {gj})")
                system.add_system(gravity[gi], (gravity[gj],))
    #print("")

    model_time = 0|units.Myr
    ax = plot_cluster(particles, model_time)

    #dt = 0.25|units.Myr
    dt = o.dt
    system.timestep = 0.5*dt

    #t_end = 10 | units.Myr
    t_end = o.t_end
    while model_time<t_end:
        model_time += dt
        system.evolve_model(model_time)
        for fi in channel["from_gr"]:
            fi.copy()
        print(model_time.in_(units.Myr))
        plot_cluster(particles, model_time, ax)
    system.stop()
    plot_cluster(particles, model_time, ax, o.figname)


