#### INITIALIZE_SSTARS
####
#### Initialize S-star orbits from observed quantities to Cartesian coordinates
####
#### Adopts the observed S-star orbits for the 27 stars observed by
#### 2009ApJ...692.1075G.
#### Observed quantities are converted to SI, and subsequently the stars are
#### put in a Keplerian orbit around the supermassive-black hole.
####
#### Dec 2022, Simon Portegies Zwart
####

from amuse.lab import *

S_name = ["S1","S2","S4","S5","S6","S8","S9","S12","S13","S14","S17","S18","S19","S21","S24","S27","S29","S31","S33","S38","S66","S67","S71","S83","S87","S96","S97"]
S_a_arcsec = [0.508,0.123,0.298,0.250,0.436,0.411,0.293,0.308,0.297,0.256,0.311,0.265,0.798,0.213,1.060,0.454,0.397,0.298,0.410,0.139,1.210,1.095,1.061,2.785,1.260,1.545,2.186]
S_ecc = [0.496,0.880,0.406,0.842,0.886,0.824,0.825,0.900,0.490,0.963,0.364,0.759,0.844,0.784,0.933,0.952,0.916,0.934,0.731,0.802,0.178,0.368,0.844,0.657,0.423,0.131,0.302]
S_inc = [120.8,135.3,77.83,143.7,86.4,74.0,81.0,31.6,25.5,99.4,96.4,116.0,73.58,54.8,106.3,92.9,122,153.8,42.9,166,135.4,139.9,76.3,123.8,142.7,126.8,114.6]
S_omra = [341.61,225.39,258.11,109,83.46,315.90,147.58,240.4,73.1,227.74,188.06,215.2,342.9,252.7,4.2,191.90,157.2,103,82.9,286,96.8,106.0,34.6,73.6,109.9,115.78,107.72]
S_Omega = [115.3,63.56,316.4,236.3,129.5,345.2,225.2,308.8,248.2,339.0,319.5,151.7,153.3,182.6,291.5,308.2,343.3,314,328.1,203,106,215.2,331.4,197.2,41.5,231.0,38]
S_tperi = [0.95,2.32,-25.6,-16.4,63,-16.2,-12.2,-4.37,4.90,0.07,-8.0,-4.0,5.1,28.1,24.9,59.7,21,13.8,-32.1,3.0,-218,-305,-354,61,-353,-376,175] | units.yr
S_Period = [132,15.8,59.5,45.7,105,96.1,58,62.5,59.2,47.3,63.2,50,260,35.8,398,112,91,59.4,96,18.9,486,419,399,1700,516,701,1180]| units.yr

def initialize_sstars(time, name, a_arcsec, ecc, inc, omra, Omega, tperi, Period, reverse_time=False):

    Rgc = 8.178 | units.kpc
    BH = Particles(1)
    BH.mass = 4.154e+6 | units.MSun
    BH.name = "SMBH"
    BH.position = (0,0,0) |units.AU 
    BH.velocity = (0,0,0) | (units.AU / units.yr)

    stars = Particles(len(name))
    stars.name = "BH    " # reserve space to store the stellar names
    for si in range(len(name)):
        S = stars[si]
        S.name = name[si]
        S.mass = 20.0 | units.MSun
        S.radius = 0 |units.RSun
        a = Rgc * a_arcsec[si] * (1 |units.AU)/(1 | units.parsec)
        print("a=", a.in_(units.au))
        from ic.orbital_elements_to_cartesian import orbital_elements_to_pos_and_vel
        S.position, S.velocity = orbital_elements_to_pos_and_vel(time, a, ecc[si], inc[si], omra[si], Omega[si], tperi[si], Period[si], BH[0].mass, S.mass)
    if reverse_time:
        for i in range(len(Sstars)):
            stars[i].velocity = -1 * stars[i].velocity
        for i in range(len(black_hole)):
            BH[i].velocity = -1 * BH[i].velocity

    return BH, stars

def plot_stars_in_XY(black_hole, SStars):
    from matplotlib import pyplot, rc
    fig = pyplot.figure(figsize=(12,12))
    font = {'size' : 20}
    rc('font', **font)
    pyplot.figaspect(1.0)
    pyplot.title("S-stars on the 1st of January 2001")
    pyplot.xlabel("X [AU]")
    pyplot.ylabel("Y [AU]")
#    fig.set.aspect(1.0)
    pyplot.figaspect(1.0)
    positions = SStars.position
    x, y, z = positions.x.value_in(units.AU), positions.y.value_in(units.AU), positions.z.value_in(units.AU)
    vx, vy, vz = SStars.vx.value_in(units.kms), SStars.vy.value_in(units.kms), SStars.vz.value_in(units.kms)
    xvx = x + vx
    yvy = y + vy
    zvz = z + vz

    pyplot.scatter([0.0],[0.0], 155000, c='w')
    pyplot.scatter([0.0],[0.0], 15500, c='w')
    pyplot.scatter([0.0],[0.0], 1550, c='w')
    pyplot.scatter([0.0], [0.0], 155, c='k')
    pyplot.scatter(x, y, 30, c='r')
    pyplot.scatter([0.0], [0.0], 30, c='k')
    pyplot.scatter([0.0], [0.0], 30, c='k')
    for i in range(len(x)):
        pyplot.arrow(x[i], y[i], vx[i], vy[i], head_width=0.05, head_length=0.1, fc='k', ec='k')
    pyplot.text(8100, 7500, "1000AU")
    pyplot.text(2000, 3450, "100AU")
    pyplot.arrow(-13000, -28000, 1000, 0, head_width=0.05, head_length=0.1, fc='k', ec='k')
    pyplot.text(-13000, -27500, "1000km/s")

#    fig.savefig("hydro_experiment_i"+str(index)+".png")
    pyplot.show()

def new_sstar_particle_set(date):
    BH, SStars = initialize_sstars(date, S_name, S_a_arcsec, S_ecc, S_inc, S_omra, S_Omega, S_tperi, S_Period)
    return BH, SStars
    
def generate_SStar_initial_conditions(time = 2001|units.yr,
                                      plot=False):
    BH, SStars = initialize_sstars(time, S_name, S_a_arcsec, S_ecc, S_inc, S_omra, S_Omega, S_tperi, S_Period)
    for si in SStars:
        print(si.name, "&", si.mass, "&", round(si.x.value_in(units.AU), 4), "&", round(si.y.value_in(units.AU), 4), "&", round(si.x.value_in(units.AU), 4), "&", round(si.vx.value_in(units.kms), 4), "&", round(si.vy.value_in(units.kms), 4), "&", round(si.vz.value_in(units.kms), 4), "\\\\")

    if plot:
        plot_stars_in_XY(BH, SStars)
        
    particles = Particles()
    BH.radius = 2*constants.G*BH.mass/constants.c**2
    SStars.radius = (SStars.mass/(1|units.MSun))*0.78 | units.RSun
    particles.add_particles(BH)
    particles.add_particles(SStars)
    return particles


def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-t", unit=units.yr,
                      dest="time", type="float", default = 2001|units.yr,
                      help="time since BC [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    particles = generate_SStar_initial_conditions(**o.__dict__)
    write_set_to_file(particles, "SStars_initial_conditions.amuse")

