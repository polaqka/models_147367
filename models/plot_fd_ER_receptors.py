
###

import numpy
import matplotlib.pyplot as plt

###

res = numpy.loadtxt(f'results/ip3r_fd_results_ER_receptors_10.dat')
ca_concs = numpy.loadtxt(f'results/ip3r_fd_ca_concs_ER_receptors_10.dat')
states =  numpy.array(['A00','A01','A10','A11','Pa','Pb','Pc','Sa','Sb','o1','o2','o3','Ia','Ib']) 
    
plt.figure(figsize = (20,10))
plt.suptitle('Receptors state after 2000 s, Fraiman & Dawson model \n [IP_3] = 10 nM, cyt. [Ca²⁺] = 70 nM', fontsize = 22)
for i in range(ca_concs.size):
    plt.subplot(4,2,i+1)
    plt.plot(states, res[i,:], 'o', markersize = 8)
    plt.ylim((-0.25,6.5))
    plt.xlabel('State', fontsize = 10)
    plt.ylabel('Mean number of receptors', fontsize = 10)
    plt.title(f'ER lumen [Ca²⁺] = {ca_concs[i]*1e9} nM', fontsize = 14)
    plt.grid(alpha = 0.3)
plt.savefig('charts/P0_ER_fd_receptors_10.jpg', format = 'jpg')
plt.close()