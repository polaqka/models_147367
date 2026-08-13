#
# Katri Hituri
#
# Script to simulate open probability of IP3R
# Model of Dawson et al. 2003
#
#
# for different cytosolic IP3 concentrations 
# and in fixed Ca²⁺ concentration [Ca²⁺] = 0.25 uM
# in steady-state


####

import ip3r_model_dli as model

import steps.rng as srng
import steps.solver as ssolver
import numpy

####



#Concentrations for IP3 in cytosol
ip_concs = numpy.array([0.001e-6, 0.01e-6, 0.05e-6,
                        0.07e-6, 0.10e-6, 0.15e-6, 0.20e-6, 0.25e-6, 
                        0.30e-6, 0.33e-6, 0.35e-6, 0.36e-6, 0.43e-6, 0.50e-6, 
                        1.00e-6, 1.25e-6, 1.50e-6, 2.00e-6, 2.50e-6, 
                        3e-6, 5.00e-6, 10.00e-6, 30e-6, 100e-6, 300e-6, 1000e-6]) # mol/l

# Solver settings
r = srng.create('mt19937', 1000)
r.initialize(2605)
sim = ssolver.Wmdirect(model.mdl, model.cell, r)

# Number of iterations (defines how many times the model is simulated)
NITER = 5000


tpnt = numpy.arange(0.0, 50.01, 0.01)
# array for simulation results
res = numpy.zeros([ip_concs.size, 2])

print('Simulating the IP3R model of Dawson et al. 2003.')

for i in range(ip_concs.size):
        print(('Round', i+1, '/', ip_concs.size))
        temp_res = numpy.zeros([NITER, tpnt.size])
        for j in range(NITER): 
                sim.reset()
                sim.setPatchSpecCount('ER_memb', 'R', 1)  # number of naive receptor
                sim.setCompSpecConc('cyt', 'IP3', ip_concs[i])
                sim.setCompSpecConc('cyt', 'Ca', 0.25e-6)  # [Ca²+] = 0.25 uM
                sim.setCompSpecClamped('cyt', 'Ca', 1)
                sim.setCompSpecClamped('cyt', 'IP3', 1)
                for t in range(tpnt.size):
                        sim.run(tpnt[t])
                        o1 = sim.getPatchSpecCount('ER_memb', 'O1')
                        o2 = sim.getPatchSpecCount('ER_memb', 'O2')
                        temp_res[j,t] = o1 + o2
        temp = numpy.mean(temp_res[:,2501:])
        res[i,0] = numpy.mean(temp)
        res[i,1] = numpy.std(temp)
numpy.savetxt('results/ip3r_dli_res_ip.dat', res)
numpy.savetxt('results/ip3r_dli_ip_concs.dat', ip_concs)

print('sim_dli_ip done')

