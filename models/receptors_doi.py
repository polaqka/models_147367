
#
# Script to simulate open probability of IP3R kinetics
# The model of Doi et al. 2005
#
# to analyze receptor state transitions in time
# under fixed [Ca²⁺] and varying [IP3]

###

import ip3r_model_doi as model

import steps.rng as srng
import steps.solver as ssolver
import numpy

####

# Solver initialization

r = srng.create('mt19937', 1000)

r.initialize(26058)

sim = ssolver.Wmdirect(model.mdl, model.cell, r)


# Number of iterations (defines how many times the model is simulated)

NITER = 100 #1500


# timepoint array

tpnt = numpy.arange(0.0, 14400.1, 0.1)


# Ca2+ concentrations in cytosol
'''
ca_concs = numpy.array([0.001e-6, 0.003e-6, 0.007e-6, 0.01e-6, 0.013e-6, 0.03e-6, 0.10e-6, 0.13e-6,
						 0.20e-6, 0.27e-6, 0.28e-6, 0.30e-6, 0.33e-6, 0.4e-6, 0.50e-6, 0.6e-6, 0.7e-6,
						 0.8e-6, 1.00e-6, 1.50e-6, 3.00e-6, 10.00e-6, 30.00e-6, 100.00e-6]) # mol/l
'''
ca_concs = numpy.array([100e-9])


ip_concs = numpy.array([10e-9, 20e-9, 30e-9, 40e-9, 100e-9, 150e-9, 
						200e-9, 500e-9,700e-9, 1e-6, 1.5e-6, 2e-6])
#ip_concs = numpy.array([10e-9])


# array for simulation results
res = numpy.zeros([ca_concs.size, tpnt.size, ip_concs.size, 7])
 


print('Simulating the IP3R model of Doi et al. 2005.')
print('You can abort the simulation by pressing Ctrl + C')

for ip in range(ip_concs.size):
	for i in range(ca_concs.size):

		print('Ca = ', ca_concs[i],', ip3 = ', ip_concs[ip])
		temp_res = numpy.zeros([NITER, tpnt.size, 7])  # temporary storage for results

		for j in range(NITER): 
	
			sim.reset()
			sim.setPatchSpecCount('ER_memb', 'R', 10)  # number of naive receptor
			sim.setCompSpecConc('cyt', 'IP3', ip_concs[ip])
			sim.setCompSpecConc('cyt', 'Ca', ca_concs[i])
			sim.setCompSpecClamped('cyt', 'Ca', 1) # Ca in cytosol is constant
			sim.setCompSpecClamped('cyt', 'IP3', 1) # IP3 in cytosol is constant
			#
			for t in range(tpnt.size):
				sim.run(tpnt[t]) # run the simulation
				R = sim.getPatchSpecCount('ER_memb','R')
				R_ip3 = sim.getPatchSpecCount('ER_memb', 'RIP3')
				R_open = sim.getPatchSpecCount('ER_memb', 'Ropen')
				RCa = sim.getPatchSpecCount('ER_memb', 'RCa')
				RCa2 = sim.getPatchSpecCount('ER_memb', 'RCa2')
				RCa3 = sim.getPatchSpecCount('ER_memb', 'RCa3')
				RCa4 = sim.getPatchSpecCount('ER_memb', 'RCa4')
				temp_res[j,t, :] = numpy.array([R, R_open, R_ip3, RCa, RCa2, RCa3, RCa4])
		# calculate the mean and standard deviation of the simulation results      
		res[i,:,ip,:] =  numpy.mean(temp_res, axis = 0 )

# save the results (means and stds)
print(res.shape)
res_reshaped = res.reshape(res.shape[0], -1)
numpy.savetxt('ip3r_doi_results_receptors.dat', res_reshaped)
numpy.savetxt('ip3r_doi_ca_receptors.dat', ca_concs)
numpy.savetxt('ip3r_doi_ip_receptors.dat', ip_concs)
