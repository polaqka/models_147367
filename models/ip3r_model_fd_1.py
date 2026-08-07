# ip3r_model.py
#
# Katri Hituri
# IP3R model SUSTAINABLE 
###########

# Import packages

import steps.model as smodel
import steps.geom as sgeom



# *MODEL*
#

mdl = smodel.Model()
volsys = smodel.Volsys('vsys', mdl)  # Volume system
surfsys = smodel.Surfsys('ssys', mdl)  # Surface system


# CHEMICAL SPECIES 
Ca = smodel.Spec('Ca', mdl)  # Calcium
IP3 = smodel.Spec('IP3', mdl)  # IP3

# IP3 receptor states
A00 = smodel.Spec('A00', mdl)  # naive state
A01 = smodel.Spec('A01', mdl)
A10 = smodel.Spec('A10', mdl)  # open state
A11 = smodel.Spec('A11', mdl)
Pa = smodel.Spec('Pa', mdl)
Pb = smodel.Spec('Pb', mdl)
Pc = smodel.Spec('Pc', mdl)
Sa =  smodel.Spec('Sa', mdl)
Sb =  smodel.Spec('Sb', mdl)
Oa =  smodel.Spec('Oa', mdl)
Ob = smodel.Spec('Ob', mdl)
Oc = smodel.Spec('Oc', mdl)
Ia = smodel.Spec('Ia', mdl)
Ib = smodel.Spec('Ib', mdl)

# REACTIONS


#1 Sa + Ca <=> Sb
reac1_f = smodel.SReac('Sa_Sb', surfsys, olhs=[Ca], slhs=[Sa], srhs=[Sb])
reac1_b = smodel.SReac('Sb_Sa', surfsys, slhs=[Sb], orhs=[Ca], srhs=[Sa])

#2 Pc <=> Sa
reac2_f = smodel.SReac('Pc_Sa', surfsys, slhs=[Pc], srhs=[Sa])
reac2_b = smodel.SReac('Sa_Pc', surfsys, slhs=[Sa], srhs=[Pc])

#3 Pb + Ca_ER <=> Pc
reac3_f = smodel.SReac('Pb_Pc', surfsys, ilhs=[Ca], slhs=[Pb], srhs=[Pc])
reac3_b = smodel.SReac('Pc_Pb', surfsys, slhs=[Pc], irhs=[Ca], srhs=[Pb])

#4 Pa <=> Pb
reac4_f = smodel.SReac('Pa_Pb', surfsys, slhs=[Pa], srhs=[Pb])
reac4_b = smodel.SReac('Pb_Pa', surfsys, slhs=[Pb], srhs=[Pa])

#5 A01 <=> Pa
reac5_f = smodel.SReac('A01_Pa', surfsys, slhs=[A01], srhs=[Pa])
reac5_b = smodel.SReac('Pa_A01', surfsys, slhs=[Pa], srhs=[A01])

#6 A00 + Ca <=> A01
reac6_f = smodel.SReac('A00_A01', surfsys, olhs=[Ca], slhs=[A00], srhs=[A01])
reac6_b = smodel.SReac('A01_A00', surfsys, slhs=[A01], orhs=[Ca], srhs=[A00]) 

#7 A00 + IP3 <=> A10
reac7_f = smodel.SReac('A00_A10', surfsys, olhs=[IP3], slhs=[A00], srhs=[A10])
reac7_b = smodel.SReac('A10_A00', surfsys, slhs=[A10], orhs=[IP3], srhs=[A00])

#8 A01 + IP3 <=> A11
reac8_f = smodel.SReac('A01_A11', surfsys, olhs=[IP3], slhs=[A01], srhs=[A11])
reac8_b = smodel.SReac('A11_A01', surfsys, slhs=[A11], orhs=[IP3], srhs=[A01])

#9 A10 + Ca <=> A11
reac9_f = smodel.SReac('A10_A11', surfsys, olhs=[Ca], slhs=[A10], srhs=[A11])
reac9_b = smodel.SReac('A11_A10', surfsys, slhs=[A11], orhs=[Ca], srhs=[A10])

#10 A11 <=> Oa
reac10_f = smodel.SReac('A11_Oa', surfsys, slhs=[A11], srhs=[Oa])
reac10_b = smodel.SReac('Oa_A11', surfsys, slhs=[Oa], srhs=[A11])

#11 Oa <=> Ob
reac11_f = smodel.SReac('Oa_Ob', surfsys, slhs=[Oa], srhs=[Ob])
reac11_b = smodel.SReac('Ob_Oa', surfsys, slhs=[Ob], srhs=[Oa])

#12 Ob + Ca_ER <=> Oc
reac12_f = smodel.SReac('Ob_Oc', surfsys, ilhs=[Ca], slhs=[Ob], srhs=[Oc])
reac12_b = smodel.SReac('Oc_Ob', surfsys, slhs=[Oc], irhs=[Ca], srhs=[Ob])

#13 Oc <=> Ia
reac13_f = smodel.SReac('Oc_Ia', surfsys, slhs=[Oc], srhs=[Ia])
reac13_b = smodel.SReac('Ia_Oc', surfsys, slhs=[Ia], srhs=[Oc])

#14 Ia + Ca <=> Ib
reac14_f = smodel.SReac('Ia_Ib', surfsys, olhs=[Ca], slhs=[Ia], srhs=[Ib])
reac14_b = smodel.SReac('Ib_Ia', surfsys, slhs=[Ib], orhs=[Ca], srhs=[Ia])

# Rate constants and their units
reac1_f.kcst = 5000e6       # 1/(Ms)
reac1_b.kcst = 20           # 1/s
reac2_f.kcst = 3000         # 1/(Ms)
reac2_b.kcst = 250          # 1/s
reac3_f.kcst = 5000e6       # 1/(Ms)
reac3_b.kcst = 150          # 1/s
reac4_f.kcst = 500          # 1/s
reac4_b.kcst = 100          # 1/s
reac5_f.kcst = 0.3          # 1/s
reac5_b.kcst = 700          # 1/s
reac6_f.kcst = 5000e6       # 1/(Ms)
reac6_b.kcst = 1            # 1/s
reac7_f.kcst = 6670e6       # 1/(Ms)
reac7_b.kcst = 0.0117       # 1/s
reac8_f.kcst = 1540e6       # 1/(Ms)
reac8_b.kcst = 18           # 1/s
reac9_f.kcst = 500e6        # 1/(Ms)
reac9_b.kcst = 667          # 1/s
reac10_f.kcst = 1800        # 1/s
reac10_b.kcst = 330         # 1/s
reac11_f.kcst = 133         # 1/s
reac11_b.kcst = 1500        # 1/s
reac12_f.kcst = 70e6        # 1/(Ms)
reac12_b.kcst = 2000        # 1/s
reac13_f.kcst = 630         # 1/s
reac13_b.kcst = 400         # 1/s
reac14_f.kcst = 60e6        # 1/(Ms)
reac14_b.kcst = 16          # 1/s

    

# GEOMETRY

cell = sgeom.Geom()

# Create the cytosol compartment
cyt = sgeom.Comp('cyt', cell)					
cyt.addVolsys('vsys')
cyt.setVol(0.1e-18)

# Create the endoplasmic reticulum lumen compartment
ER_lumen = sgeom.Comp('ER_lumen', cell)					
ER_lumen.addVolsys('vsys')
ER_lumen.setVol(0.1e-18)

# Create the er membrane (ER_lumen inner, cyt outer)
ER_memb = sgeom.Patch('ER_memb', cell, ER_lumen, cyt)		
ER_memb.addSurfsys('ssys')
ER_memb.setArea(0.21544369e-12)



