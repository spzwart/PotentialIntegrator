import numpy as np
from amuse.datamodel import Particles
from amuse.units import units, quantities
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

from potentials.point_particle_potential import PointParticlePotential
from potentials.plummer_potential import PlummerPotential

class CompositeGravityCode(object):
    def __init__(self, converter, potential=PointParticlePotential()):
        self.model_time = 0 | units.Myr
        self.particles = Particles()
        self.converter = converter
        self.potential = potential

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
        m = self.particles[0].mass
        x = self.particles[0].x
        y = self.particles[0].y
        z = self.particles[0].z
        eps = self.potential.epsilon2
        return m*self.potential.get_potential_at_point(eps,x,y,z)[0]
    
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

def get_kinetic_energy(gravity):
    Ek = 0 | units.erg
    for gi in gravity:
        Ek += gi.kinetic_energy
    return Ek
        
def get_potential_energy(gravity):
    Ep = 0 | units.erg
    for gi in gravity:
        Ep += gi.potential_energy
    return Ep
    
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
        #potential = PointParticlePotential()
        potential = PlummerPotential()
        gravity.append(CompositeGravityCode(converter, potential))    
        gravity[-1].add_particle(particles[pi].as_set())
        channel["to_gr"].append(particles.new_channel_to(gravity[-1].particles))
        channel["from_gr"].append(gravity[-1].particles.new_channel_to(particles))

        
    system=bridge.Bridge(verbose=False)
    for gi in range(len(gravity)):
        for gj in range(len(gravity)):
            if gi != gj:
                print(f"add_system({gi}, ({gj},))")
                system.add_system(gravity[gi], (gravity[gj],))
    #print("")

    model_time = 0|units.Myr
    ax = plot_cluster(particles, model_time)

    dt = o.dt
    system.timestep = 0.25*dt

    Ek0 = get_kinetic_energy(gravity)
    Ep0 = get_potential_energy(gravity)
    
    t_end = o.t_end
    while model_time<t_end:
        model_time += dt
        system.evolve_model(model_time)
        for fi in channel["from_gr"]:
            fi.copy()
        Ek = get_kinetic_energy(gravity)
        Ep = get_potential_energy(gravity)
        dE = (Ek-Ek0) + (Ep-Ep0)
        print(model_time.in_(units.Myr), dE/(Ek+Ep))
        plot_cluster(particles, model_time, ax)
        
    system.stop()
    plot_cluster(particles, model_time, ax, o.figname)


