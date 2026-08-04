#
# Katri Hituri
#
# Script to simulate open probability of IP3R
# The model of Doi et al. 2005
#

####

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
NITER = 1500

# timepoint array
tpnt = numpy.arange(0.0, 40.01, 0.01)

#Concentrations for Ca in cytosol
ca_concs = numpy.array([0.01e-6, 0.02e-6, 0.05e-6, 0.07e-6, 0.10e-6, 0.15e-6, 0.20e-6, 0.25e-6, 0.28e-6, 0.30e-6, 0.33e-6, 0.35e-6, 0.36e-6, 0.38e-6, 0.43e-6, 0.50e-6, 1.00e-6, 1.50e-6, 2.00e-6, 2.50e-6, 5.00e-6]) # mol/l

# array for simulation results
res = numpy.zeros([ca_concs.size, 2])



print 'Simulating the IP3R model of Doi et al. 2005.'
print 'You can abort the simulation by pressing Ctrl + C'

for i in xrange(ca_concs.size):

	print 'Round', i+1, '/', ca_concs.size
	temp_res = numpy.zeros([NITER, tpnt.size]) # temporary storage for results

	for j in xrange(NITER): 
   
	        sim.reset()
	        sim.setPatchCount('ER_memb', 'R', 1) # number of naive receptor
	        sim.setCompConc('cyt', 'IP3', 10e-6) # [IP3] = 10 uM
	        sim.setCompConc('cyt', 'Ca', ca_concs[i])
	        sim.setCompClamped('cyt', 'Ca', 1) # Ca in cytosol is constant
	        sim.setCompClamped('cyt', 'IP3', 1) # IP3 in cytosol is constant
	    
	        for t in xrange(tpnt.size):
	            sim.run(tpnt[t]) # run the simulation
	            temp_res[j,t] = sim.getPatchCount('ER_memb', 'Ropen')
        
	# calculate the mean and standard deviation of the simulation results      
    	temp = numpy.mean(temp_res[:,2501:]) # take only into account results after 25 s
    	res[i,0] = numpy.mean(temp, 0)
    	res[i,1] = numpy.std(temp, 0)


# save the results (means and stds)
numpy.savetxt('ip3r_doi_op_res.dat', res)
numpy.savetxt('ip3r_doi_op_ca_concs.dat', ca_concs)

print res[:,0]




