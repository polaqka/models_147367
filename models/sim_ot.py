#
# Katri Hituri
#
# Script to simulate open probability of IP3R
# Model of Othmer and Tang 1993
#
# for different cytosolic Ca²⁺ concentrations 
# and in different IP3 concentrations
# in steady-state


####

import ip3r_model_ot as model

import steps.rng as srng
import steps.solver as ssolver
import numpy

####


#Concentrations for Ca and IP3 in cytosol
ca_concs = numpy.array([0.001e-6, 0.003e-6, 0.007e-6, 0.01e-6, 0.013e-6, 0.03e-6, 
						0.05e-6, 0.06e-6, 0.07e-6, 0.08e-6, 0.09e-6, 0.11e-6,
						0.13e-6, 0.15e-6, 0.20e-6, 0.30e-6, 0.50e-6, 0.8e-6, 1.00e-6,
						 1.50e-6, 3.00e-6, 10.00e-6]) # mol/l
#ip_concs = numpy.array([0.02e-6, 0.1e-6,0.2e-6,1e-6,2e-6,10e-6])  # original
ip_concs = numpy.array([10e-9, 20e-9, 30e-9, 40e-9, 100e-9, 150e-9, 200e-9])  # smaller


# Solver settings
r = srng.create('mt19937', 1000)
r.initialize(26058)
sim = ssolver.Wmdirect(model.mdl, model.cell, r)

NITER = 12000

tpnt = numpy.arange(0.0, 50.01, 0.01)

# array for simulation results
res = numpy.zeros([ca_concs.size, ip_concs.size, 2])

print('Simulating the IP3R model of Othmer and Tang 1993.')

for ip in range(ip_concs.size):
	for i in range(ca_concs.size):
		#print('Round', i+1, '/', ca_concs.size)
		temp_res = numpy.zeros([NITER, tpnt.size])
		for j in range(NITER):
			sim.reset()
			sim.setPatchSpecCount('ER_memb', 'R000', 1) 
			sim.setCompSpecConc('cyt', 'IP3', ip_concs[ip])
			sim.setCompSpecClamped('cyt', 'IP3', 1)
			sim.setCompSpecConc('cyt', 'Ca', ca_concs[i])
			sim.setCompSpecClamped('cyt', 'Ca', 1)
			#
			for t in range(tpnt.size):
				sim.run(tpnt[t])
				temp_res[j,t] = sim.getPatchSpecCount('ER_memb', 'Ropen')
		
		temp = numpy.mean(temp_res[:,2501:], axis = 1)
		res[i,ip,0] = numpy.mean(temp)
		res[i,ip,1] = numpy.std(temp)

res_reshaped = res.reshape(res.shape[0], -1)
numpy.savetxt('results/ip3r_ot_results.dat', res_reshaped)
numpy.savetxt('results/ip3r_ot_ca_c.dat', ca_concs)
numpy.savetxt('results/ip3r_ot_ip_c.dat', ip_concs)

print('sim_ot done')