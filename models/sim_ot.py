#
# Katri Hituri
#
# Script to simulate open probability of IP3R
# Model of Othmer and Tang 1993




####
import ip3r_model_ot as model

import steps.rng as srng
import steps.solver as ssolver
import numpy

####



#Concentrations for Ca in cytosol
ca_concs = numpy.array([0.001e-6, 0.003e-6, 0.007e-6, 0.01e-6, 0.013e-6, 0.03e-6, 0.05e-6, 0.06e-6, 0.07e-6, 0.08e-6, 0.09e-6, 0.11e-6, 0.13e-6, 0.15e-6, 0.20e-6, 0.30e-6, 0.50e-6, 0.8e-6, 1.00e-6, 1.50e-6, 3.00e-6, 10.00e-6]) # mol/l

# Solver settings
r = srng.create('mt19937', 1000)
r.initialize(26058)
sim = ssolver.Wmdirect(model.mdl, model.cell, r)

NITER = 1  # 12000

tpnt = numpy.arange(0.0, 50.01, 0.01)

# array for simulation results
res = numpy.zeros([ca_concs.size, tpnt.size])

print('Simulating the IP3R model of Othmer and Tang 1993.')

for i in range(ca_concs.size):
	print('Round', i+1, '/', ca_concs.size)
	temp_res = numpy.zeros([NITER, tpnt.size])
	for j in range(NITER):
		sim.reset()
		sim.setPatchSpecCount('ER_memb', 'R000', 1) 
		sim.setCompSpecConc('cyt', 'IP3', 10e-6) # [IP3] = 10 uM
		sim.setCompSpecClamped('cyt', 'IP3', 1)
		sim.setCompSpecConc('cyt', 'Ca', ca_concs[i])
		sim.setCompSpecClamped('cyt', 'Ca', 1)
		#
		for t in range(tpnt.size):
			sim.run(tpnt[t])
			temp_res[j,t] = sim.getPatchSpecCount('ER_memb', 'Ropen')
	
	temp = numpy.mean(temp_res[:,2501:])
	res[i,0] = numpy.mean(temp)
	res[i,1] = numpy.std(temp)


numpy.savetxt('ip3r_ot_op_res.dat', res)
numpy.savetxt('ip3r_ot_op_ca_concs.dat', ca_concs)

print(res[:,0])


