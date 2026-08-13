#
# Pola Kukielka
#
# Plot open probability of IP3R
# as a function of cytosolic Ca²⁺ concentration 
# in fixed IP_3 concentration [IP_3] = 10 uM
# for all 4 models


###

import numpy
import matplotlib.pyplot as plt

###

models = ['ot', 'dli', 'fd', 'doi']
colors = ['g', 'b', 'r', 'm']

for i, model in enumerate(models):
    if model == 'fd':
        res = numpy.loadtxt(f'results/ip3r_{model}_res_ca_1.dat')
        ca_concs = numpy.loadtxt(f'results/ip3r_{model}_ca_concs_1.dat')
    else:
	    res = numpy.loadtxt(f'results/ip3r_{model}_res_ca.dat')
	    ca_concs = numpy.loadtxt(f'results/ip3r_{model}_ca_concs.dat')
        
    plt.plot(ca_concs*1e6, res[:,0],  linewidth = 2, color = colors[i], label = f'model {model}')
    plt.xscale('log')
    plt.axis([1e-3,1e2,0,0.8])
    plt.xlabel('[Ca²⁺] uM', fontsize = 10)
    plt.ylabel('P_o', fontsize = 10)
    plt.suptitle('Open probability of IP_3R as a function of cytosolic [Ca²⁺]', fontsize = 14)
    plt.title('[IP_3] = 10 uM', fontsize = 12)
    plt.grid(alpha = 0.3)
    #plt.legend(loc = 'upper left')

plt.savefig('charts/P0_Ca_all.jpg', format = 'jpg')
plt.close()