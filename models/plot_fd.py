###
import numpy as np
import matplotlib.pyplot as plt

op_fd = np.loadtxt('ip3r_fd_op_res.dat')
op_fd_ca_concs = np.loadtxt('ip3r_fd_op_ca_concs.dat')

op_fd_1 = np.loadtxt('ip3r_fd_op_res_1.dat')
op_fd_ca_concs_1 = np.loadtxt('ip3r_fd_op_ca_concs_1.dat')

plt.figure()
plt.plot(op_fd_ca_concs*1e6, op_fd, linewidth = 2, color = 'r', label = 'original')
plt.plot(op_fd_ca_concs_1*1e6, op_fd_1, linewidth = 2, color = 'g', label = 'sustainable')
plt.xscale('log')
plt.axis([1e-3,1e2,0,0.5])
plt.xlabel('[Ca^{2+}] uM', fontsize = 16)
plt.ylabel('P_o', fontsize = 16)
plt.suptitle('Open probability of IP_3R as a function of cytosolic [Ca^{2+}]')
plt.title('Fraiman and Dawson model')
plt.grid(alpha = 0.3)
#plt.legend(loc = 'upper left')

plt.savefig('P0_Ca.jpg', format = 'jpg')

plt.close()