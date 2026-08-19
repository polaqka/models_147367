
###

import numpy
import matplotlib.pyplot as plt

###

res = numpy.loadtxt(f'results/ip3r_fd_res_ER.dat')
ca_concs = numpy.loadtxt(f'results/ip3r_fd_ca_concs_ER.dat')
     
plt.plot(ca_concs*1e9, res[:,0],  linewidth = 2)
#plt.xscale('log')
plt.axis([-10,260,0.1,0.18])
plt.xlabel('ER lumen [Ca²⁺] nM', fontsize = 10)
plt.ylabel('P_o', fontsize = 10)
plt.suptitle('Open probability of IP_3R, Fraiman & Dawson', fontsize = 14)
plt.title('[IP_3] = 10 nM, cyt. [Ca²⁺] = 70 nM', fontsize = 10)
plt.grid(alpha = 0.3)
#plt.legend(loc = 'upper left')

plt.savefig('charts/P0_ER_fd.jpg', format = 'jpg')
plt.close()