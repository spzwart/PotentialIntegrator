#PotentialIntegrator

python make_initial_conditions.py -I SStars
python make_initial_conditions.py -I Binary

python gravity_potential.py -f binary.amuse
 PotentialIntegrator"


Timings for N=10 Rvir=1pc Plummer for 10Myr with 0.1Myr steps
code	    	   time
gravity_pure 	   1.3s
gravity_bridge	   1.8s
gravity_potential  15.2s

Timings for two 10MSun stars in a 10pc circular orbit
code	    	   time
gravity_pure 	   1.3s
gravity_bridge	   1.9
gravity_potential  1.5s