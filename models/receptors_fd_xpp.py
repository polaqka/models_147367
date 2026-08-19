
#
# Script to simulate open probability of IP3R kinetics
# The model of Fraiman and Dawson 2004
#
# to analyze receptor state transitions in time
# under fixed [Ca²⁺] and [IP3] 
# in order to compare with differential equations

###

import ip3r_model_fd_1 as model

import steps.rng as srng
import steps.solver as ssolver
import numpy

####


# Ca2+ and IP3 concentrations in cytosol
ca_concs = numpy.array([70e-9]) # [Ca_cyt] = 70 nM
ip_concs = numpy.array([10e-9]) # [ip3] = 10 nM

# Solver settings
r = srng.create('mt19937', 1000)
r.initialize(2605)
sim = ssolver.Wmdirect(model.mdl, model.cell, r)

# Number of iterations (defines how many times the model is simulated)
NITER = 1000  # 750

dt = 0.0001
tpnt = numpy.arange(0.0, 0.020 + dt, dt) # 20 ms	

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
			sim.setCompSpecConc('ER_lumen', 'Ca', 100e-9) # Ca_lumen = 100 nM
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
		print(f'po 10ms: {res[i,ip,100,:]}')
		print(f'po 20ms: {res[i,ip,200,:]}')

# save the results (means and stds)
print(res.shape)
res_reshaped = res.reshape(res.shape[0], -1)
numpy.savetxt('results/ip3r_fd_results_receptors_xpp.dat', res_reshaped)
numpy.savetxt('results/ip3r_fd_ca_receptors_xpp.dat', ca_concs)
numpy.savetxt('results/ip3r_fd_ip_receptors_xpp.dat', ip_concs)
