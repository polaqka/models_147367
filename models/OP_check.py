#
# De Young and Keizer model (1992) 
# compare open probability of IP3R 
# from the equation and steps simulations


import numpy
import matplotlib.pyplot as plt

loaded_res = numpy.loadtxt(f'results/ip3r_dyk_results.dat')
ca_concs = numpy.loadtxt(f'results/ip3r_dyk_ca_c.dat')
ip_concs = numpy.loadtxt(f'results/ip3r_dyk_ip_c.dat')

res = loaded_res.reshape((ca_concs.size, ip_concs.size, 2))


# EQUATION 
# rate constants
a1 = 4e8       # 1/(Ms)
b1 = 52        # 1/s
a2 = 2e5       # 1/(Ms)
b2 = 0.21      # 1/s
a3 = 4e8       # 1/(Ms)
b3 = 377.36    # 1/s
a4 = 2e5       # 1/Ms
b4 = 0.0289    # 1/s
a5 = 2e7       # 1/(Ms)
b5 = 1.6468    # 1/s

d1 = 0.13e-6   # b1/a1
d2 = 1.049e-6  # b2/a2
d3 = 934.4e-9  # b3/a3
d4 = 144.5e-9  # b4/a4
d5 = 82.34e-9  # b5/a5

'''
Kd1 = 145e-9  # M
Kd2 = 542e-9  # M
IP3 = 15e-9 # M
d1 = Kd1 - IP3
d3 = (Kd2 - IP3)*(1+d2)- d1*d2
'''

P_o = numpy.zeros((ca_concs.size, ip_concs.size))

for ip in range(ip_concs.size):
    for ca in range(ca_concs.size):
        P_o[ca,ip] = ((ca_concs[ca]*ip_concs[ip]*d2)
                        /((ca_concs[ca]*ip_concs[ip]+ip_concs[ip]*d2 + d1*d2 
                        + ca_concs[ca]*d3)*(ca_concs[ca]+d5)))**3


plt.figure(figsize = (30,12))
plt.suptitle('Open probability of IP_3R as a function of [Ca²⁺] for different [IP_3] \n \n De Young & Keizer model', fontsize = 22)
plt.subplot(1,2,1)

xmax = numpy.max(ca_concs*1e6)
xmin = numpy.min(ca_concs*1e6)
#ymax = numpy.max(res[:,:,0])


for ip in range(ip_concs.size):
    plt.plot(ca_concs*1e6, P_o[:,ip], label = f'ip3 = {ip_concs[ip]*1e6} uM')
plt.xscale('log')
ymax = numpy.max(P_o)
plt.axis([xmin,xmax,0,1.3*ymax])
plt.xlabel('cyt. [Ca²⁺] uM', fontsize = 16)
plt.ylabel('P_o', fontsize = 16)
plt.title('from equation', fontsize = 18)
plt.grid(alpha = 0.3)
plt.legend(loc = 'upper left')


plt.subplot(1,2,2)
for ip in range(ip_concs.size):
    plt.plot(ca_concs*1e6, res[:,ip,0], label = f'ip3 = {ip_concs[ip]*1e6} uM')
plt.xscale('log')
ymax = numpy.max(res[:,:,0])
plt.axis([xmin,xmax,0,1.3*ymax])
plt.xlabel('cyt. [Ca²⁺] uM', fontsize = 16)
plt.ylabel('P_o', fontsize = 16)
plt.title('from STEPS simulation', fontsize = 18)
plt.grid(alpha = 0.3)
plt.legend(loc = 'upper left')

plt.savefig('charts/P0_dyk_compare.jpg', format = 'jpg')
plt.close()