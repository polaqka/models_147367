

import model_dyk as model

import steps.rng as srng
import steps.solver as ssolver
import numpy

ca_concs = numpy.array([200e-9])  # , 0.1e-6, 0.01e-6])

ip_concs = numpy.array([2e-6])  # ([10e-9, 20e-9, 30e-9, 40e-9, 100e-9, 150e-9, 200e-9])

r = srng.create('mt19937', 1000)
r.initialize(2605)
sim = ssolver.Wmdirect(model.mdl, model.cell, r)

NITER = 750
dt = 0.001
tpnt = numpy.arange(0.0, 30 + dt, dt)

mean_open_time = numpy.zeros([ca_concs.size, ip_concs.size, 2])

print('Simulating mean open time for IP3R model (De Young and Keizer 1992)...')

for ip in range(ip_concs.size):
    for i in range(ca_concs.size):

        all_open_times = []

        for j in range(NITER):
            sim.reset()
            sim.setPatchSpecCount('ER_memb', 'S000', 1)
            sim.setCompSpecConc('cyt', 'IP3', ip_concs[ip])
            sim.setCompSpecClamped('cyt', 'IP3', 1)
            sim.setCompSpecConc('cyt', 'Ca', ca_concs[i])
            sim.setCompSpecClamped('cyt', 'Ca', 1)
            #sim.setCompSpecConc('ER_lumen', 'Ca', 150e-6)
            #sim.setCompSpecClamped('ER_lumen', 'Ca', 1)

            # open state = 1, closed state = 0
            is_open = numpy.zeros(tpnt.size)

            for t in range(tpnt.size):
                sim.run(tpnt[t])
                S100 = sim.getPatchSpecCount('ER_memb', 'S110')
                if S100 > 0: is_open[t] = 1

            diff = numpy.diff(is_open)   #[int(20.0 / dt):])
            starts = numpy.where(diff == 1)[0] + 1
            ends = numpy.where(diff == -1)[0] + 1

            if len(starts) > 0 and len(ends) > 0:
                if ends[0] < starts[0]:
                    ends = ends[1:]
                if len(starts) > len(ends):
                    starts = starts[:len(ends)]

                open_times = (ends - starts) * dt
                all_open_times.extend(open_times)


        if len(all_open_times) > 0:
            mean_open_time[i, ip, 0] = numpy.mean(all_open_times)
            mean_open_time[i, ip, 1] = numpy.std(all_open_times)
        else:
            mean_open_time[i, ip, 0] = 0
            mean_open_time[i, ip, 1] = 0

        print(f'Done Ca {i+1}/{ca_concs.size}, IP3 {ip+1}/{ip_concs.size}')

res_reshaped = mean_open_time.reshape(mean_open_time.shape[0], -1)
numpy.savetxt('results/ip3r_dyk_mean_open_time.dat', res_reshaped)
numpy.savetxt('results/ip3r_dyk_open_times.dat', all_open_times)
