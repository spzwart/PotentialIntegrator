#PotentialIntegrator

This is a simple N-body integrator using potentials rather than
particles.  The potentials are coupled together to form an N-potential
system (rather than an N-body system).  The potentials are
subsequently bridged to form a self-consistent system than can be
integrated. At the moment, integration is realized using second-order
symplectic Verlet, but there is no particular reason why this couldn't
be a higher even order time-summetric symplectic integrator.


The following directories:
.         # this directory
ic        # repository for initial conditions
plot      # plotting routine(s)
potential # potentials

Directory: "./"
 README.md			  # This file
 gravity_pure.py		  # Performs 4th order N-body integration
 gravity_potential.py		  # The actual potential integrator
 read_parameters.py		  # parameter reader
 make_initial_conditions.py  	  # initial conditions generator

Directory: "./potential"
 point_particle_potential.py      # Point-particle potential (mimics N-body)
 plummer_potential.py             # Plummer potential


Directory: "./plot"
 plot_cluster.py		  # plot simulation result
 
Directory: "./ic"
 initialize_sstars.py		  # generate initial conditions for S-stars
 orbital_elements_to_cartesian.py # helper routine to convert coordinates


Example of how to use the code
# generate initial conditions
python make_initial_conditions.py -I SStars # to generate the S-star cluster
python make_initial_conditions.py -I Binary # to generate a simple binary

# run the code as pure N-body (4th order)
python gravity_pure.py -f binary.amuse -t 5.e+7 --dt 1.e+6

# Or run the code with potentials instead of particles (2nd order)
python gravity_potential.py -f binary.amuse -t 5.e+7 --dt 1.e+6

For running the S-star cluster for 10
python make_initial_conditions.py -I SStars 
python gravity_potential.py -f SStars.amuse -t 10 dt 0.1 &

Check the result with the pure N-body code:
python gravity_pure.py -f SStars.amuse -t 10 dt 0.1 &

For running a 10 particle Plummer the sphere
python make_initial_conditions.py -I Plummer -N 10 --Rvir 1
python gravity_potential.py -f plummer.amuse -t 1.e+7 --dt 1.e+5

To check the same initial setup with the direct N-body run:
python gravity_pure.py -f plummer.amuse -t 10 --dt 0.05


Timings for a few runs.  The direct (pure) N-body code is considerably
faster for any number of particles. This is not a surprise. It lacks
uptimization at the moment. Performance could be improved by adopting
the fast-kick routines in AMUSE, but that would require GPU hardware
to be available.

Timings for N=10 Rvir=1pc Plummer for 10Myr with 0.1Myr steps
code	    	   time
gravity_pure 	   1.3s
gravity_potential  15.2s

Timings for two 10MSun stars in a 10pc circular orbit
code	    	   time
gravity_pure 	   1.3s
gravity_potential  1.5s