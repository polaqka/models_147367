#
# Pola Kukielka
#
# Plot open probability of IP3R
# as a function of cytosolic [IP3] concentration 
# in fixed Ca²⁼ concentration [Ca²⁼] = 0.25 uM
# for Fraiman and Dawson model
# original and sustainable k_rates (microscopic reversibility)


###

import numpy as np
import matplotlib.pyplot as plt

###


fd = np.loadtxt('results/ip3r_fd_res_ip.dat')
fd_ip_concs = np.loadtxt('results/ip3r_fd_ip_concs.dat')

fd_1 = np.loadtxt('results/ip3r_fd_res_ip_1.dat')
fd_ip_concs_1 = np.loadtxt('results/ip3r_fd_ip_concs_1.dat')

plt.figure()
plt.plot(fd_ip_concs*1e6, fd[:,0], linewidth = 2, color = 'r', label = 'original')
plt.plot(fd_ip_concs_1*1e6, fd_1[:,0], linewidth = 2, color = 'g', label = 'sustainable')
plt.xscale('log')
plt.axis([1e-3,1e3,0,0.5])
plt.xlabel('cytosolic [IP3] uM', fontsize = 10)
plt.ylabel('P_o', fontsize = 10)
plt.suptitle('Open probability of IP_3R, Fraiman & Dawson model', fontsize = 14)
plt.title('microscopic reversibility comparison', fontsize = 10)
plt.grid(alpha = 0.3)
plt.legend(loc = 'upper right')

plt.savefig('charts/P0_IP3_fd.jpg', format = 'jpg')

plt.close()