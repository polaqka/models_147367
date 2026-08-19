#
# Pola Kukielka
#
# Plot open probability of IP3R
# as a function of cytosolic Ca²⁺ concentration 
# in fixed IP3 concentration [IP3] = 10 uM
# for Doi et al. model
# original and extended


###

import numpy as np
import matplotlib.pyplot as plt

###


doi = np.loadtxt('results/ip3r_doi_res_ip.dat')
doi_ip_concs = np.loadtxt('results/ip3r_doi_ip_concs.dat')

doi_2 = np.loadtxt('results/ip3r_doi_2_res_ip.dat')
doi_ip_concs_2 = np.loadtxt('results/ip3r_doi_2_ip_concs.dat')

plt.figure()
plt.plot(doi_ip_concs*1e6, doi[:,0], linewidth = 2, label = 'original')
plt.plot(doi_ip_concs_2*1e6, doi_2[:,0], linewidth = 2, label = 'extended')
plt.xscale('log')
plt.axis([1e-2,1e3,0,0.5])
plt.xlabel('[IP3] uM', fontsize = 10)
plt.ylabel('P_o', fontsize = 10)
plt.suptitle('Open probability of IP_3R, Doi et al. model', fontsize = 14)
plt.title('comparison', fontsize = 10)
plt.grid(alpha = 0.3)
plt.legend(loc = 'upper right')

plt.savefig('charts/P0_IP3_doi.jpg', format = 'jpg')

plt.close()