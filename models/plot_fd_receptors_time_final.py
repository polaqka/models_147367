
###

import numpy
import matplotlib.pyplot as plt

###

res = numpy.loadtxt(f'results/ip3r_fd_results_receptors_xpp_final.dat')
ca_concs = numpy.loadtxt(f'results/ip3r_fd_ca_concs_receptors_xpp_final.dat')
states =  numpy.array(['A00','A01','A10','A11','Pa','Pb','Pc','Sa','Sb','o1','o2','o3','Ia','Ib']) 
    
plt.figure(figsize = (15,10))
plt.suptitle('Receptors state after 20 ms, Fraiman & Dawson model \n [IP_3] = 10 nM, cyt. [Ca²⁺_cyt] = 70 nM, [Ca_lum] = 100 nM', fontsize = 22)
for i in range(ca_concs.size):
    plt.plot(states, res[:], 'o', markersize = 8)
    plt.ylim((-0.25,6.5))
    plt.xlabel('State', fontsize = 12)
    plt.ylabel('Mean number of receptors', fontsize = 12)
    plt.grid(alpha = 0.3)
plt.savefig('charts/P0_fd_receptors_xpp_final.jpg', format = 'jpg')
plt.close()