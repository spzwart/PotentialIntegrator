
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

if __name__ == "__main__":
    o, arguments = new_option_parser().parse_args()

    particles = read_set_from_file(o.filename, close_file=True)
    m = particles.mass.sum()
    r = particles.position.length().sum()/len(particles)
    converter=nbody_system.nbody_to_si(m, r)

    gravity = Ph4(converter)
    gravity.particles.add_particles(particles)
    gravity.parameters.timestep_parameter = 0.01
    c2gry=particles.new_channel_to(gravity.particles)
    c2fr=gravity.particles.new_channel_to(particles)

    Ek0 = gravity.kinetic_energy
    Ep0 = gravity.potential_energy
    
    #ax = plot_cluster(potentials)
    model_time = 0|units.Myr
    ax = plot_cluster(particles, model_time)

    dt = o.dt
    t_end = o.t_end
    while model_time<t_end:
        model_time += dt
        gravity.evolve_model(model_time)
        c2fr.copy()
        Ek = gravity.kinetic_energy
        Ep = gravity.potential_energy
        dE = (Ek-Ek0) + (Ep-Ep0)
        print(model_time.in_(units.Myr), dE/(Ek+Ep))
        plot_cluster(particles, model_time, ax)
    gravity.stop()
    plot_cluster(particles, model_time, ax, o.figname)

