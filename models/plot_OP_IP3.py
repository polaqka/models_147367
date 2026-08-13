#
# Pola Kukielka
#
# Plot open probability of IP3R
# as a function of cytosolic IP_3 concentration 
# in fixed Ca²⁺concentration [Ca²⁺] = 0.25 uM
# for all 4 models


###

import numpy
import matplotlib.pyplot as plt

###

models = ['ot', 'dli', 'fd', 'doi']
colors = ['g', 'b', 'r', 'm']

for i, model in enumerate(models):
    res = numpy.loadtxt(f'results/ip3r_{model}_res_ip.dat')
    ip_concs = numpy.loadtxt(f'results/ip3r_{model}_ip_concs.dat')
        
    plt.plot(ip_concs*1e6, res[:,0],  linewidth = 2, color = colors[i], label = f'model {model}')
    plt.xscale('log')
    plt.axis([1e-3,1e3,0,0.5])
    plt.xlabel('[Ca²⁺] uM', fontsize = 10)
    plt.ylabel('P_o', fontsize = 10)
    plt.suptitle('Open probability of IP_3R as a function of cytosolic [IP_3]', fontsize = 14)
    plt.title('[Ca²⁺] = 0.25 uM', fontsize = 12)
    plt.grid(alpha = 0.3)
    #plt.legend(loc = 'upper left')

plt.savefig('charts/P0_IP3_all.jpg', format = 'jpg')
plt.close()