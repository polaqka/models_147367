
#
# Script to simulate open probability of IP3R
# The model of Fraiman and Dawson 2004
#
# SUSTAINABLE
#
# for different ER_lumen Ca²⁺ concentrations 
# in fixed cytosolic Ca²⁺ and IP3 concentrations 
# in steady-state


####

import ip3r_model_fd_1 as model

import steps.rng as srng
import steps.solver as ssolver
import numpy

####


# Ca2+ concentrations in ER
ca_concs = numpy.array([0, 1e-9,10e-9, 50e-9, 100e-9, 150e-9,200e-9, 250e-9]) # mol/l

# Solver settings
r = srng.create('mt19937', 1000)
r.initialize(2605)
sim = ssolver.Wmdirect(model.mdl, model.cell, r)

# Number of iterations (defines how many times the model is simulated)
NITER = 750

dt = 0.01
tpnt = numpy.arange(0.0, 30+dt, dt)

# array for simulation results
res = numpy.zeros([ca_concs.size, 2])

print('Simulating the IP3R model of Fraiman and Dawson 2004.')

for i in range(ca_concs.size):

	#print('Round', i+1, '/', ca_concs.size)
	temp_res = numpy.zeros([NITER, tpnt.size]) # temporary storage for results

	for j in range(NITER): 
		sim.reset()
		sim.setPatchSpecCount('ER_memb', 'A00', 1) # number of naive receptor
		sim.setCompSpecConc('cyt', 'IP3', 10e-9) # [IP3] = 10 nM
		sim.setCompSpecClamped('cyt', 'IP3', 1)
		sim.setCompSpecConc('cyt', 'Ca', 70e-9)  # [Ca²⁺] = 70 nm
		sim.setCompSpecClamped('cyt', 'Ca', 1)
		sim.setCompSpecConc('ER_lumen', 'Ca', ca_concs[i]) 
		sim.setCompSpecClamped('ER_lumen', 'Ca', 1)
		#
		for t in range(tpnt.size):
			sim.run(tpnt[t])
			o1 = sim.getPatchSpecCount('ER_memb', 'Oa')
			o2 = sim.getPatchSpecCount('ER_memb', 'Ob')
			o3 = sim.getPatchSpecCount('ER_memb', 'Oc')  
			temp_res[j,t] = o1 + o2 + o3 

	# calculate the mean and standard deviation of the simulation results   
	temp = numpy.mean(temp_res[:, 2001:], axis = 1)
	res[i,0] = numpy.mean(temp)
	res[i,1] = numpy.std(temp)


# save the results (means and stds)
numpy.savetxt('results/ip3r_fd_res_ER.dat', res)
numpy.savetxt('results/ip3r_fd_ca_concs_ER.dat', ca_concs)

print('sim_fd_ER done')

