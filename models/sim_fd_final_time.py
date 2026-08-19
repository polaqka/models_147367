
#
# Script to simulate open probability of IP3R
# The model of Fraiman and Dawson 2004
#
# SUSTAINABLE
#
# for different ER_lumen Ca²⁺ concentrations 
# in fixed cytosolic Ca²⁺ and IP3 concentrations 
# in time


####

import ip3r_model_fd_1 as model

import steps.rng as srng
import steps.solver as ssolver
import numpy

####


# Ca2+ concentrations in ER
ca_concs = numpy.array([100e-9])

# Solver settings
r = srng.create('mt19937', 1000)
r.initialize(2605)
sim = ssolver.Wmdirect(model.mdl, model.cell, r)

# Number of iterations (defines how many times the model is simulated)
NITER = 100 #750

dt = 0.001
tpnt = numpy.arange(0.0, 0.20+dt, dt)

# array for simulation results
res = numpy.zeros([ca_concs.size, 14])

print('Simulating the IP3R model of Fraiman and Dawson 2004.')

for i in range(ca_concs.size):

	print('Round', i+1, '/', ca_concs.size)
	temp_res = numpy.zeros([NITER, 14]) # temporary storage for results

	for j in range(NITER): 
		sim.reset()
		sim.setPatchSpecCount('ER_memb', 'A00', 10) # number of naive receptor
		sim.setCompSpecConc('cyt', 'IP3', 10e-9) # [IP3] = 10 nM
		sim.setCompSpecClamped('cyt', 'IP3', 1)
		sim.setCompSpecConc('cyt', 'Ca', 70e-9)  # [Ca²⁺] = 70 nM
		sim.setCompSpecClamped('cyt', 'Ca', 1)
		sim.setCompSpecConc('ER_lumen', 'Ca', ca_concs[i]) # [Ca_lum] = 100 nM
		sim.setCompSpecClamped('ER_lumen', 'Ca', 1)
		#
		for t in range(tpnt.size):
			sim.run(tpnt[t])
			A00 = sim.getPatchSpecCount('ER_memb', 'A00')
			A01 = sim.getPatchSpecCount('ER_memb', 'A01')
			A10 = sim.getPatchSpecCount('ER_memb', 'A10')
			A11 = sim.getPatchSpecCount('ER_memb', 'A11')
			Pa = sim.getPatchSpecCount('ER_memb', 'Pa')
			Pb = sim.getPatchSpecCount('ER_memb', 'Pb')
			Pc = sim.getPatchSpecCount('ER_memb', 'Pc')
			Sa = sim.getPatchSpecCount('ER_memb', 'Sa')
			Sb = sim.getPatchSpecCount('ER_memb', 'Sb')  				
			o1 = sim.getPatchSpecCount('ER_memb', 'Oa')
			o2 = sim.getPatchSpecCount('ER_memb', 'Ob')
			o3 = sim.getPatchSpecCount('ER_memb', 'Oc') 
			Ia = sim.getPatchSpecCount('ER_memb', 'Ia')
			Ib = sim.getPatchSpecCount('ER_memb', 'Ib')
        # receptors state after 20 ms
		temp_res[j,:] =  numpy.array([A00,A01,A10,A11,Pa,Pb,Pc,Sa,Sb,o1,o2,o3,Ia,Ib])
		#print(temp_res[j,:])
		#print(numpy.sum(temp_res[j,:]))

    # calculate the mean from the iterations   
	res[i,:] =  numpy.mean(temp_res, axis = 0)
	#print(res[i,:])
	#print(numpy.sum(res[i,:]))

# save the results (means and stds)
numpy.savetxt('results/ip3r_fd_results_receptors_xpp_final.dat', res)
numpy.savetxt('results/ip3r_fd_ca_concs_receptors_xpp_final.dat', ca_concs)

print('sim_fd_xpp_time done')

