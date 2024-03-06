from amuse.lab import *
import numpy as np
from matplotlib import pyplot as plt

def plot_cluster(particles,
                 model_time=0|units.Myr,
                 ax1=None,
                 figfile=""):

    if ax1==None:
        plt.rcParams.update({'font.size': 20})
        fig, ax1 = plt.subplots(1, 1, figsize=(10, 10))
        ax1.yaxis.tick_right()
        ax1.yaxis.tick_left()
        ax1.xaxis.tick_top()
        ax1.xaxis.tick_bottom()
        ax1.minorticks_on() # switch on the minor ticks
        ax1.locator_params(nbins=3)
        ax1.set_xlim(-10, 10)
        ax1.set_ylim(-10, 10)
        ax1.axis("equal")

    #print(f"t={potentials[0].model_time.in_(units.Myr)}, x={x}")
    s = 1 + model_time.value_in(units.Myr)
    from matplotlib import cm
    c = cm.rainbow(np.linspace(0, 1, len(particles)))
    ax1.scatter(particles.x.value_in(units.pc),
                particles.y.value_in(units.pc),
                s=s, c=c)
    #for i in range(len(m)):
    #    plot_circle(x[i], y[i], r[i], ax1)

    if len(figfile)>0:
        plt.savefig(figfile)
        plt.show()
        
    return ax1

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("--fig", 
                      dest="filename", 
                      default = "binary.amuse",
                      help="initial conditions datafile [%default]")
    result.add_option("-f", 
                      dest="figname", 
                      default = "figure.pdf",
                      help="figure filename [%default]")
    return result


if __name__ == "__main__":
    o, arguments = new_option_parser().parse_args()

    particles = read_set_from_file(o.filename, close_file=True)
    plot_cluster(particles,
                 figfile=o.figname)

