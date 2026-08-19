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


doi = np.loadtxt('results/ip3r_doi_res_ca.dat')
doi_ca_concs = np.loadtxt('results/ip3r_doi_ca_concs.dat')

doi_2 = np.loadtxt('results/ip3r_doi_2_res_ca.dat')
doi_ca_concs_2 = np.loadtxt('results/ip3r_doi_2_ca_concs.dat')

plt.figure()
plt.plot(doi_ca_concs*1e6, doi[:,0], linewidth = 2, label = 'original')
plt.plot(doi_ca_concs_2*1e6, doi_2[:,0], linewidth = 2, label = 'extended')
plt.xscale('log')
plt.axis([0.5e-2,1e1,0,0.5])
plt.xlabel('cytosolic [Ca²⁺] uM', fontsize = 10)
plt.ylabel('P_o', fontsize = 10)
plt.suptitle('Open probability of IP_3R, Doi et al. model', fontsize = 14)
plt.title('comparison', fontsize = 10)
plt.grid(alpha = 0.3)
plt.legend(loc = 'upper right')

plt.savefig('charts/P0_Ca_doi.jpg', format = 'jpg')

plt.close()