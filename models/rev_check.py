###
# microscopic reversibility check for two models

# fraiman and dawson

reac1_f = 5000e6       # 1/(Ms)
reac1_b = 20           # 1/s
reac2_f = 3000         # 1/(Ms)
reac2_b = 250          # 1/s
reac3_f = 5000e6       # 1/(Ms)
reac3_b = 150          # 1/s
reac4_f = 500          # 1/s
reac4_b = 100          # 1/s
reac5_f = 0.3          # 1/s
reac5_b = 700          # 1/s
reac6_f = 5000e6       # 1/(Ms)
reac6_b = 1            # 1/s
reac7_f = 6670e6       # 1/(Ms)
reac7_b = 200          # 1/s
reac8_f = 1540e6       # 1/(Ms)
reac8_b = 18           # 1/s
reac9_f = 500e6        # 1/(Ms)
reac9_b = 667          # 1/s
reac10_f = 1800        # 1/s
reac10_b = 330         # 1/s
reac11_f = 133         # 1/s
reac11_b = 1500        # 1/s
reac12_f = 70e6        # 1/(Ms)
reac12_b = 2000        # 1/s
reac13_f = 630         # 1/s
reac13_b = 400         # 1/s
reac14_f = 60e6        # 1/(Ms)
reac14_b = 16          # 1/s

print('Fraiman and Dawson')
print(f'Reversibility check: {reac6_b * reac8_b/reac6_f/reac8_f == 
                                reac7_b * reac9_b/reac7_f/reac9_f}')
print(f'better k_b_7 = {reac6_b * reac8_b/reac6_f/reac8_f/reac9_b*reac7_f*reac9_f} 1/s')


# Dawson at al.

reac1_f = 1           # 1/s
reac1_b = 100         # 1/s

reac2_f = 4000e6      # 1/(Ms) 
reac2_b = 1000        # 1/s

reac3_f = 3000e6      # 1/(Ms)
reac3_b = 2000        # 1/s

reac4_f = 2000e6      # 1/Ms
reac4_b = 3000        # 1/s

reac5_f = 1000e6      # 1/Ms
reac5_b = 4000        # 1/s

reac6_f = 400e6       # 1/Ms
reac6_b = 10          # 1/s

reac7_f = 300e6       # 1/Ms
reac7_b = 20          # 1/s

reac8_f = 200e6       # 1/Ms
reac8_b = 30          # 1/s

reac9_f = 100e6       # 1/Ms
reac9_b = 40          # 1/s

reac10_f = 1          # 1/s
reac10_b = 10         # 1/s

reac11_f = 1          # 1/s
reac11_b = 1          # 1/s

reac12_f = 10         # 1/s
reac12_b = 1          # 1/s

reac13_f = 10         # 1/s
reac13_b = 0.1        # 1/s

reac15_f = 100e6      # 1/Ms
reac15_b = 10         # 1/s

reac18_f = 1e6        # 1/Ms
reac18_b = 0.1        # 1/s

reac19_f = 10e6       # 1/Ms
reac19_b = 0.01       # 1/s

    

index = [[1,6,2,10], [10,7,3,11], [11,8,4,12],[12,9,5,13]]

print(' ')
print('Dawson')
print(f'Reversibility check 1: {reac1_b * reac6_b/reac1_f/reac6_f == 
                                reac2_b * reac10_b/reac2_f/reac10_f}')
print(f'Reversibility check 2: {reac10_b * reac7_b/reac10_f/reac7_f == 
                                reac3_b * reac11_b/reac3_f/reac11_f}')
print(f'Reversibility check 3: {reac11_b * reac8_b/reac11_f/reac8_f == 
                                reac4_b * reac12_b/reac4_f/reac12_f}')
print(f'Reversibility check 4: {reac12_b * reac9_b/reac12_f/reac9_f == 
                                reac5_b * reac13_b/reac5_f/reac13_f}')
print(f'Reversibility check 5: {reac1_b * reac19_b/reac1_f/reac19_f == 
                                reac18_b/reac18_f}')

print(f'better k_b_19 = {reac1_b * reac18_f/reac1_f/reac18_b*reac19_f}')
