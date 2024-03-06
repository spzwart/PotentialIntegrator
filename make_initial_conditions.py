import numpy as np
from amuse.lab import *
from amuse.ic.plummer import new_plummer_model
from amuse.ic.salpeter import new_salpeter_mass_distribution

from ic.initialize_sstars import generate_SStar_initial_conditions

def new_PlummerModel(N, Rvir):
    converter=nbody_system.nbody_to_si(N|units.MSun, Rvir)
    particles = new_plummer_model(N, converter)
    particles.mass = new_salpeter_mass_distribution(len(particles))
    particles.scale_to_standard(converter)
    particles.move_to_center()
    return particles

def new_binary_star():
    N = 2
    mass = 10 | units.MSun
    radius = 0.001|units.pc
    converter=nbody_system.nbody_to_si(mass.sum(), radius)
    particles = Particles(2)
    vc = np.sqrt(constants.G*2*mass/(1|units.pc))
    particles[0].position = [1,0,0] | units.pc
    particles[0].velocity = [0,1,0] *vc
    particles[1].position = [0,0,0] | units.pc
    particles[1].velocity = [0,0,0] *vc
    particles.mass = mass
    particles.radius = radius
    particles.move_to_center()
    return particles

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-I", 
                      dest="conditions", 
                      default = "Binary",
                      help="initial conditions routine [%default]")
    result.add_option("-N", 
                      dest="N",  type = "int",
                      default = 10,
                      help="Number of stars [%default]")
    result.add_option("-R", unit=units.pc,
                      dest="Rvir",  type = "float",
                      default = 1|units.pc,
                      help="Cluster virial radius [%default]")
    result.add_option("-f", 
                      dest="filename", 
                      default = "",
                      help="initial conditions datafile [%default]")
    return result

if __name__ == "__main__":
    o, arguments = new_option_parser().parse_args()
    set_printing_strategy("custom", 
                          preferred_units = [units.MSun,
                                             units.au,
                                             units.Myr], 
                          precision = 6, prefix = "", 
                          separator = " [", suffix = "]")


    if "Plummer" in o.conditions:
        particles = new_PlummerModel(o.N, o.Rvir)
        filename = "plummer.amuse"
    elif "Binary" in o.conditions:
        particles = new_binary_star()
        filename = "binary.amuse"
    elif "SStars":
        particles = generate_SStar_initial_conditions()
        filename = "SStars.amuse"
    else:
        print("No initial conditions specified.")
        print("Stop.")
        exit(-1)

    if len(o.filename)>0:
        filename = o.filename
        
    write_set_to_file(particles, filename)
