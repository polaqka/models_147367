#
# Katri Hituri
#
# Script to simulate open probability of IP3R
# The model of Fraiman and Dawson 2004
#
# 
# for different cytosolic IP3 concentrations 
# and in fixed Ca²⁺ concentration [Ca²⁺] = 0.25 uM
# in steady-state



####

import ip3r_model_fd as model

import steps.rng as srng
import steps.solver as ssolver
import numpy

####


# IP3 concentrations in cytosol
ip_concs = numpy.array([0.001e-6, 0.003e-6, 0.007e-6, 0.01e-6, 0.013e-6, 0.03e-6, 0.10e-6, 0.13e-6,
						 0.20e-6, 0.27e-6, 0.28e-6, 0.30e-6, 0.33e-6, 0.4e-6, 0.50e-6, 0.6e-6, 0.7e-6,
						 0.8e-6, 1.00e-6, 1.50e-6, 3.00e-6, 10.00e-6, 30.00e-6, 100.00e-6, 500e-6, 1000e-6]) # mol/l

# Solver settings
r = srng.create('mt19937', 1000)
r.initialize(2605)
sim = ssolver.Wmdirect(model.mdl, model.cell, r)

# Number of iterations (defines how many times the model is simulated)
NITER = 750

tpnt = numpy.arange(0.0, 30.01, 0.01)

# array for simulation results
res = numpy.zeros([ip_concs.size, tpnt.size])

print('Simulating the IP3R model of Fraiman and Dawson 2004.')

for i in range(ip_concs.size):

	print('Round', i+1, '/', ip_concs.size)
	temp_res = numpy.zeros([NITER, tpnt.size]) # temporary storage for results

	for j in range(NITER): 
		sim.reset()
		sim.setPatchSpecCount('ER_memb', 'A00', 1) # number of naive receptor
		sim.setCompSpecConc('cyt', 'IP3', ip_concs[i])
		sim.setCompSpecClamped('cyt', 'IP3', 1)
		sim.setCompSpecConc('cyt', 'Ca',0.25e-6)  # [Ca ²+] = 0.25 uM
		sim.setCompSpecClamped('cyt', 'Ca', 1)
		sim.setCompSpecConc('ER_lumen', 'Ca', 150e-6) 
		sim.setCompSpecClamped('cyt', 'Ca', 1)
		#
		for t in range(tpnt.size):
			sim.run(tpnt[t])
			o1 = sim.getPatchSpecCount('ER_memb', 'Oa')
			o2 = sim.getPatchSpecCount('ER_memb', 'Ob')
			o3 = sim.getPatchSpecCount('ER_memb', 'Oc')  
			temp_res[j,t] = o1 + o2 + o3 

	# calculate the mean and standard deviation of the simulation results   
	temp = numpy.mean(temp_res[:,2001:], axis = 1) # take only into account results after 20 s
	res[i,0] = numpy.mean(temp)
	res[i,1] = numpy.std(temp)


# save the results (means and stds)
numpy.savetxt('results/ip3r_fd_res_ip.dat', res)
numpy.savetxt('results/ip3r_fd_ip_concs.dat', ip_concs)


print('sim_fd_ip done')