
###

import numpy
import matplotlib.pyplot as plt

###


model_names = {'ot' : 'Othmer & Tang', 'dli': 'Dawson, Lea & Irvine',
				 'fd': 'Fraiman & Dawson', 'doi': 'Doi et al.'}
models = ['ot', 'dli', 'fd', 'doi']
colors = ['tab:blue', 'tab:red', 'tab:orange', 'tab:green']

for i, model in enumerate(models):
    open_times = numpy.loadtxt(f'results/ip3r_{model}_open_times.dat')
    plt.title('Distributions of IP_3R open times', fontsize = 12)
    plt.hist(open_times * 1e3, bins = 50, density = True, 
                histtype = 'step', label = model_names[model], color = colors[i])
    plt.xlabel('open time, ms', fontsize = 12)
    if model == 'ot' : plt.xlim((0,3500))
    if model == 'dli' : plt.xlim((0,100))
    if model == 'fd' : plt.xlim((0,30))
    if model == 'doi' : plt.xlim((0,6))
    plt.grid(alpha = 0.3)
    plt.legend(loc = 'upper right')
    print(f'{model} done')
    plt.savefig(f'charts/open_times_{model}.jpg', format = 'jpg')
    plt.close() 