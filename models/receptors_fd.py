
#
# Script to simulate open probability of IP3R kinetics
# The model of Fraiman and Dawson 2004
#
# to analyze receptor state transitions in time
# under fixed [Ca²⁺] and varying [IP3]

###

import ip3r_model_fd_1 as model

import steps.rng as srng
import steps.solver as ssolver
import numpy

####


# Ca2+ and IP3 concentrations in cytosol
'''
ca_concs = numpy.array([0.001e-6, 0.003e-6, 0.007e-6, 0.01e-6, 0.013e-6, 0.03e-6, 0.10e-6, 0.13e-6,
						 0.20e-6, 0.27e-6, 0.28e-6, 0.30e-6, 0.33e-6, 0.4e-6, 0.50e-6, 0.6e-6, 0.7e-6,
						 0.8e-6, 1.00e-6, 1.50e-6, 3.00e-6, 10.00e-6, 30.00e-6, 100.00e-6]) # mol/l
'''
ip_concs = numpy.array([10e-9, 20e-9, 30e-9, 40e-9, 100e-9, 150e-9, 
						200e-9, 500e-9,700e-9, 1e-6, 1.5e-6, 2e-6])
ca_concs = numpy.array([70e-9]) # [Ca_cyt] = 70 nM

# Solver settings
r = srng.create('mt19937', 1000)
r.initialize(2605)
sim = ssolver.Wmdirect(model.mdl, model.cell, r)

# Number of iterations (defines how many times the model is simulated)
NITER = 750

dt = 0.01
tpnt = numpy.arange(0.0, 20 + dt, dt) 	

# array for simulation results
res = numpy.zeros([ca_concs.size, ip_concs.size, tpnt.size, 14])
 
print('Simulating the IP3R model of Fraiman and Dawson 2004.')


for ip in range(ip_concs.size):
	for i in range(ca_concs.size):

		print(f'Ca = {ca_concs[i]}, ip3 = {ip_concs[ip]}')
		temp_res = numpy.zeros([NITER, tpnt.size, 14]) # temporary storage for results

		for j in range(NITER): 
			sim.reset()
			sim.setPatchSpecCount('ER_memb', 'A00', 10) # number of naive receptor
			sim.setCompSpecConc('cyt', 'IP3', ip_concs[ip])
			sim.setCompSpecClamped('cyt', 'IP3', 1)
			sim.setCompSpecConc('cyt', 'Ca', ca_concs[i])
			sim.setCompSpecClamped('cyt', 'Ca', 1)
			sim.setCompSpecConc('ER_lumen', 'Ca', 150e-6) # Ca_lumen = 150 uM
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
				temp_res[j,t,:] =  numpy.array([A00,A01,A10,A11,Pa,Pb,Pc,Sa,Sb,o1,o2,o3,Ia,Ib])
		 
		res[i,ip,:,:] =  numpy.mean(temp_res, axis = 0)
		print(f'na początku: {res[i,ip,0,:]}')

# save the results (means and stds)
print(res.shape)
res_reshaped = res.reshape(res.shape[0], -1)
numpy.savetxt('results/ip3r_fd_results_receptors_20.dat', res_reshaped)
numpy.savetxt('results/ip3r_fd_ca_receptors_20.dat', ca_concs)
numpy.savetxt('results/ip3r_fd_ip_receptors_20.dat', ip_concs)
