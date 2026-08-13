#
# Katri Hituri
#
# Script to simulate open probability of IP3R
# The model of Fraiman and Dawson 2004
#
# for different cytosolic Ca²⁺ concentrations 
# and in different IP3 concentrations
# in steady-state
# 
# SUSTAINABLE	

####

import ip3r_model_fd_1 as model

import steps.rng as srng
import steps.solver as ssolver
import numpy

####


# Ca2+ concentrations in cytosol
ca_concs = numpy.array([0.001e-6, 0.003e-6, 0.007e-6, 0.01e-6, 0.013e-6, 0.03e-6, 0.10e-6, 0.13e-6,
						 0.20e-6, 0.27e-6, 0.28e-6, 0.30e-6, 0.33e-6, 0.4e-6, 0.50e-6, 0.6e-6, 0.7e-6,
						 0.8e-6, 1.00e-6, 1.50e-6, 3.00e-6, 10.00e-6, 30.00e-6, 100.00e-6]) # mol/l

#ip_concs = numpy.array([0.02e-6, 0.1e-6,0.2e-6,1e-6,2e-6,10e-6])  # original
ip_concs = numpy.array([10e-9, 20e-9, 30e-9, 40e-9, 100e-9, 150e-9, 200e-9])  # smaller

# Solver settings
r = srng.create('mt19937', 1000)
r.initialize(2605)
sim = ssolver.Wmdirect(model.mdl, model.cell, r)

# Number of iterations (defines how many times the model is simulated)
NITER = 750

tpnt = numpy.arange(0.0, 30.01, 0.01)

# array for simulation results
res = numpy.zeros([ca_concs.size, ip_concs.size, 2])

print('Simulating the IP3R model of Fraiman and Dawson 2004.')


for ip in range(ip_concs.size):
	for i in range(ca_concs.size):

		#print('Round :', i+1, '/', ca_concs.size, ', ip3 = ', ip_concs[ip])
		temp_res = numpy.zeros([NITER, tpnt.size]) # temporary storage for results

		for j in range(NITER): 
			sim.reset()
			sim.setPatchSpecCount('ER_memb', 'A00', 1) # number of naive receptor
			sim.setCompSpecConc('cyt', 'IP3', ip_concs[ip])
			sim.setCompSpecClamped('cyt', 'IP3', 1)
			sim.setCompSpecConc('cyt', 'Ca', ca_concs[i])
			sim.setCompSpecClamped('cyt', 'Ca', 1)
			sim.setCompSpecConc('ER_lumen', 'Ca', 150e-6) 
			sim.setCompSpecClamped('ER_lumen', 'Ca', 1)
			#
			for t in range(tpnt.size):
				sim.run(tpnt[t])
				o1 = sim.getPatchSpecCount('ER_memb', 'Oa')
				o2 = sim.getPatchSpecCount('ER_memb', 'Ob')
				o3 = sim.getPatchSpecCount('ER_memb', 'Oc')  
				temp_res[j,t] = o1 + o2 + o3 

		# calculate the mean and standard deviation of the simulation results   
		temp = numpy.mean(temp_res[:,2001:], axis = 1) # take only into account results after 20 s
		res[i,ip,0] = numpy.mean(temp)
		res[i,ip,1] = numpy.std(temp)


# save the results (means and stds)
res_reshaped = res.reshape(res.shape[0], -1)
numpy.savetxt('results/ip3r_fd_results.dat', res_reshaped)
numpy.savetxt('results/ip3r_fd_ca_c.dat', ca_concs)
numpy.savetxt('results/ip3r_fd_ip_c.dat', ip_concs)

print('sim_fd done')