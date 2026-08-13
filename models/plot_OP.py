#
# Pola Kukielka
#
# Plot open probability of IP3R
# as a function of cytosolic Ca²⁺ concentration 
# in different IP_3 concentrations 
# for all 4 models


###

import numpy
import matplotlib.pyplot as plt

###

model_names = {'ot' : 'Othmer & Tang', 'dli': 'Dawson, Lea & Irvine',
				 'fd': 'Fraiman & Dawson', 'doi': 'Doi et al.'}
models = ['ot', 'dli', 'fd', 'doi']

plt.figure(figsize = (15,12))

for i, model in enumerate(models):
	loaded_res = numpy.loadtxt(f'results/ip3r_{model}_results.dat')
	ca_concs = numpy.loadtxt(f'results/ip3r_{model}_ca_c.dat')
	ip_concs = numpy.loadtxt(f'results/ip3r_{model}_ip_c.dat')

	res = loaded_res.reshape((ca_concs.size, ip_concs.size, 2))
	#print(res.shape)

	plt.subplot(2,2,i+1)
	for ip in range(ip_concs.size):
		if ip == 2: plt.plot(ca_concs*1e6, res[:,ip,0], label = f'ip3 = {int(ip_concs[ip]*1e9)+1} nM')
		else:  plt.plot(ca_concs*1e6, res[:,ip,0], label = f'ip3 = {int(ip_concs[ip]*1e9)} nM')
	plt.xscale('log')
	xmax = numpy.max(ca_concs*1e6)
	xmin = numpy.min(ca_concs*1e6)
	ymax = numpy.max(res[:,:,0])
	plt.axis([xmin,xmax,0,1.3*ymax])
	plt.xlabel('[Ca²⁺] uM', fontsize = 12)
	plt.ylabel('P_o', fontsize = 12)
	plt.suptitle('Open probability of IP_3R as a function of [Ca²⁺] for different [IP_3]', fontsize = 22)
	plt.title(model_names[model], fontsize = 16)
	plt.grid(alpha = 0.3)
	plt.legend(loc = 'upper left')
	print(f'{model} done')
plt.savefig('charts/P0.jpg', format = 'jpg')

plt.close()