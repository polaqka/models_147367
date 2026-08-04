#
# Katri Hituri
#
# Script to simulate open probability of IP3R
# The model of Fraiman and Dawson 2004
#
# 


####

import ip3r_model_fd as model

import steps.rng as srng
import steps.solver as ssolver
import numpy

####


# Ca2+ concentrations in cytosol
ca_concs = numpy.array([0.001e-6, 0.003e-6, 0.007e-6, 0.01e-6, 0.013e-6, 0.03e-6, 0.10e-6, 0.13e-6, 0.20e-6, 0.27e-6, 0.28e-6, 0.30e-6, 0.33e-6, 0.4e-6, 0.50e-6, 0.6e-6, 0.7e-6, 0.8e-6, 1.00e-6, 1.50e-6, 3.00e-6, 10.00e-6, 30.00e-6, 100.00e-6]) # mol/l

# Solver settings
r = srng.create('mt19937', 1000)
r.initialize(2605)
sim = ssolver.Wmdirect(model.mdl, model.cell, r)

# Number of iterations (defines how many times the model is simulated)
NITER = 750

tpnt = numpy.arange(0.0, 30.01, 0.01)

# array for simulation results
res = numpy.zeros([ca_concs.size, tpnt.size])

print 'Simulating the IP3R model of Fraiman and Dawson 2004.'

for i in xrange(ca_concs.size):

	print 'Round', i+1, '/', ca_concs.size
	temp_res = numpy.zeros([NITER, tpnt.size]) # temporary storage for results

	for j in xrange(NITER): 
		sim.reset()
        	sim.setPatchCount('ER_memb', 'A00', 1) # number of naive receptor
        	sim.setCompConc('cyt', 'IP3', 10e-6) # [IP3] = 10 uM
        	sim.setCompClamped('cyt', 'IP3', 1)
        	sim.setCompConc('cyt', 'Ca', ca_concs[i])
        	sim.setCompClamped('cyt', 'Ca', 1)
        	sim.setCompConc('ER_lumen', 'Ca', 150e-6) 
        	sim.setCompClamped('cyt', 'Ca', 1)

    
        	for t in xrange(tpnt.size):
            		sim.run(tpnt[t])
            		o1 = sim.getPatchCount('ER_memb', 'Oa')
            		o2 = sim.getPatchCount('ER_memb', 'Ob')
            		o3 = sim.getPatchCount('ER_memb', 'Oc')  
            		temp_res[j,t] = o1 + o2 + o3 

	# calculate the mean and standard deviation of the simulation results   
	temp = numpy.mean(temp_res[:,2001:]) # take only into account results after 20 s
	res[i,0] = numpy.mean(temp, 0)
    	res[i,1] = numpy.std(temp, 0)


# save the results (means and stds)
numpy.savetxt('ip3r_fd_op_res.dat', res)
numpy.savetxt('ip3r_fd_op_ca_concs.dat', ca_concs)

print res[:,0]

