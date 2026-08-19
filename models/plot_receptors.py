#
# Pola Kukielka
#
# Plot open probability of receptor states transitions in time
# uncomment the right one

###

import numpy
import matplotlib.pyplot as plt

###

# Doi et al.

loaded_res = numpy.loadtxt('results/ip3r_doi_results_receptors.dat')
ca_concs = numpy.loadtxt('results/ip3r_doi_ca_receptors.dat')
ip_concs = numpy.loadtxt('results/ip3r_doi_ip_receptors.dat')


res = loaded_res.reshape((1, 144001, 12, 7))  # doi
print(res.shape)

states = numpy.array(['R', 'R_open', 'R_ip3', 'RCa', 'RCa2', 'RCa3', 'RCa4'])  # doi             
t = numpy.arange(0,res.shape[1]/100, 1/100)


for s, state in enumerate(states):
    plt.figure(figsize = (30,10))
    plt.suptitle('Doi et al. - state '+ state, fontsize = 14)
    ymax = numpy.max(res[0,:,:,s])
    for ip in range(ip_concs.size):
        plt.subplot(6,2,ip+1)
        plt.plot(t,res[0,:,ip,s],  color = 'darkorange')
        plt.title('ip3 = ' + str(ip_concs[ip]))
        plt.axis([-10,1500,0,ymax+2])
        plt.grid(alpha = 0.3)
    #plt.supxlabel('Time [s]')
    #plt.supylabel('P_o')
    plt.tight_layout()
    plt.savefig(f'charts/rec_doi/state_{state}_doi.jpg', format = 'jpg')
    print(f'state {state} saved')

plt.close()



# Fraiman & Dawson


loaded_res = numpy.loadtxt('results/ip3r_fd_results_receptors_20.dat')
ca_concs = numpy.loadtxt('results/ip3r_fd_ca_receptors_20.dat')
ip_concs = numpy.loadtxt('results/ip3r_fd_ip_receptors_20.dat')

res = loaded_res.reshape((1, 2002, 12, 14)) # fd
print(res.shape)
           
states =  numpy.array(['A00','A01','A10','A11','Pa','Pb','Pc','Sa','Sb','o1','o2','o3','Ia','Ib'])  # fd
t = numpy.arange(0,res.shape[1]/100, 1/100)

for s, state in enumerate(states):
    plt.figure(figsize = (30,10))
    plt.suptitle('Fraiman & Dawson - state '+ state + ' po 20s', fontsize = 14)
    ymax = numpy.max(res[0,:,:,s])
    for ip in range(ip_concs.size):
        plt.subplot(6,2,ip+1)
        plt.plot(t,res[0,:,ip,s])
        plt.title('ip3 = ' + str(ip_concs[ip]))
        plt.axis([-1,21,0,ymax])
        plt.grid(alpha = 0.3)
    #plt.supxlabel('Time [s]')
    #plt.supylabel('P_o')
    plt.tight_layout()
    plt.savefig(f'charts/rec_20/state_{state}_fd_20.jpg', format = 'jpg')
    print(f'state {state} saved')

plt.close()


# Fraiman & Dawson - xpp


loaded_res = numpy.loadtxt('results/ip3r_fd_results_receptors_xpp.dat')
ca_concs = numpy.loadtxt('results/ip3r_fd_ca_receptors_xpp.dat')
ip_concs = numpy.loadtxt('results/ip3r_fd_ip_receptors_xpp.dat')

res = loaded_res.reshape((1, 1, 201, 14)) # fd
print(res.shape)
           
states =  numpy.array(['A00','A01','A10','A11','Pa','Pb','Pc','Sa','Sb','o1','o2','o3','Ia','Ib'])  # fd
t = numpy.arange(0,res.shape[2]/10000, 1/10000)

plt.figure(figsize = (30,20))
plt.suptitle('Fraiman & Dawson - [ip3] = 10 nM, [Ca_cyt] = 70 nM, [Ca_lum] = 100 nM', fontsize = 14)
for s, state in enumerate(states):
    plt.subplot(7,2,s+1)
    ymax = numpy.max(res[0,0,:,s])
    plt.plot(t,res[0,0,:,s])
    plt.title('state: ' + state)
    plt.axis([-0.0001,0.0201,0,ymax])
    plt.grid(alpha = 0.3)
    plt.xlabel('Time [s]')
    plt.ylabel('P_o')
    plt.tight_layout()
    print(f'state {state} done')
plt.savefig(f'charts/states_fd_xpp.jpg', format = 'jpg')
plt.close()

