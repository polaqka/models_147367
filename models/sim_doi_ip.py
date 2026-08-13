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

#Concentrations for IP3 in cytosol
ip_concs = numpy.array([0.01e-6, 0.02e-6, 0.05e-6, 0.07e-6, 0.10e-6, 0.15e-6, 
						0.20e-6, 0.25e-6, 0.28e-6, 0.30e-6, 0.33e-6, 0.35e-6, 
						0.36e-6, 0.38e-6, 0.43e-6, 0.50e-6, 1.00e-6, 1.50e-6, 
						2.00e-6, 2.50e-6, 5.00e-6, 10e-6, 30e-6,100e-6, 500e-6, 1000e-6]) # mol/l

# array for simulation results
res = numpy.zeros([ip_concs.size, 2])



print('Simulating the IP3R model of Doi et al. 2005.')
print('You can abort the simulation by pressing Ctrl + C')

for i in range(ip_concs.size):

	print('Round', i+1, '/', ip_concs.size)
	temp_res = numpy.zeros([NITER, tpnt.size])  # temporary storage for results

	for j in range(NITER): 
   
		sim.reset()
		sim.setPatchSpecCount('ER_memb', 'R', 1)  # number of naive receptor
		sim.setCompSpecConc('cyt', 'IP3', ip_concs[i]) 
		sim.setCompSpecConc('cyt', 'Ca', 0.25e-6)  # [Ca²+] = 0.25 uM
		sim.setCompSpecClamped('cyt', 'Ca', 1) # Ca in cytosol is constant
		sim.setCompSpecClamped('cyt', 'IP3', 1) # IP3 in cytosol is constant
		#
		for t in range(tpnt.size):
			sim.run(tpnt[t]) # run the simulation
			temp_res[j,t] = sim.getPatchSpecCount('ER_memb', 'Ropen')
        
	# calculate the mean and standard deviation of the simulation results      
	temp = numpy.mean(temp_res[:,2501:]) # take only into account results after 25 s
	res[i,0] = numpy.mean(temp)
	res[i,1] = numpy.std(temp)


# save the results (means and stds)
numpy.savetxt('ip3r_doi_res_ip.dat', res)
numpy.savetxt('ip3r_doi_ip_concs.dat', ip_concs)

print(res[:,0])




