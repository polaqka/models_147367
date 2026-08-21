

import steps.model as smodel
import steps.geom as sgeom



# *MODEL*

mdl = smodel.Model()
volsys = smodel.Volsys('vsys', mdl)  # Volume system
surfsys = smodel.Surfsys('ssys', mdl)  # Surface system

# CHEMICAL SPECIES 
Ca = smodel.Spec('Ca', mdl)  # Calcium
IP3 = smodel.Spec('IP3', mdl)  # IP3

# IP3 receptor states
S000 = smodel.Spec('S000', mdl)  # naive state
S001 = smodel.Spec('S001', mdl)
S010 = smodel.Spec('S010', mdl)
S100 = smodel.Spec('S100', mdl)
S011 = smodel.Spec('S011', mdl)
S101 = smodel.Spec('S101', mdl)
S110 = smodel.Spec('S110', mdl)  # open state
S111 = smodel.Spec('S111', mdl)

# REACTIONS


#1 S000 + IP3 <=> S100
reac1_f = smodel.SReac('S000_S100', surfsys, olhs=[IP3], slhs=[S000], srhs=[S100])
reac1_b = smodel.SReac('S100_S000', surfsys, slhs=[S100], orhs=[IP3], srhs=[S000])

#2 S100 + Ca <=> S101
reac2_f = smodel.SReac('S100_S101', surfsys, olhs=[Ca], slhs=[S100], srhs=[S101])
reac2_b = smodel.SReac('S101_S100', surfsys, slhs=[S101], orhs=[Ca], srhs=[S100])

#3 S001 + IP3 <=> S101
reac3_f = smodel.SReac('S001_S101', surfsys, olhs=[IP3], slhs=[S001], srhs=[S101])
reac3_b = smodel.SReac('S101_S001', surfsys, slhs=[S101], orhs=[IP3], srhs=[S001])

#4 S000 + Ca <=> S001
reac4_f = smodel.SReac('S000_S001', surfsys, olhs=[Ca], slhs=[S000], srhs=[S001])
reac4_b = smodel.SReac('S001_S000', surfsys, slhs=[S001], orhs=[Ca], srhs=[S000])

#5 S010 + IP3 <=> S110
reac5_f = smodel.SReac('S010_S110', surfsys, olhs=[IP3], slhs=[S010], srhs=[S110])
reac5_b = smodel.SReac('S110_S010', surfsys, slhs=[S110], orhs=[IP3], srhs=[S010])

#6 S110 + Ca <=> S111
reac6_f = smodel.SReac('S110_S111', surfsys, olhs=[Ca], slhs=[S110], srhs=[S111])
reac6_b = smodel.SReac('S111_S110', surfsys, slhs=[S111], orhs=[Ca], srhs=[S110])

#7 S011 + IP3 <=> S111
reac7_f = smodel.SReac('S011_S111', surfsys, olhs=[IP3], slhs=[S011], srhs=[S111])
reac7_b = smodel.SReac('S111_S011', surfsys, slhs=[S111], orhs=[IP3], srhs=[S011])

#8 S010 + Ca <=> S011
reac8_f = smodel.SReac('S010_S011', surfsys, olhs=[Ca], slhs=[S010], srhs=[S011])
reac8_b = smodel.SReac('S011_S010', surfsys, slhs=[S011], orhs=[Ca], srhs=[S010])

#9 S000 + Ca <=> S010
reac9_f = smodel.SReac('S000_S010', surfsys, olhs=[Ca], slhs=[S000], srhs=[S010])
reac9_b = smodel.SReac('S010_S000', surfsys, slhs=[S010], orhs=[Ca], srhs=[S000])

#10 S100 + Ca <=> S110
reac10_f = smodel.SReac('S100_S110', surfsys, olhs=[Ca], slhs=[S100], srhs=[S110])
reac10_b = smodel.SReac('S110_S100', surfsys, slhs=[S110], orhs=[Ca], srhs=[S100])

#11 S001 + Ca <=> S011
reac11_f = smodel.SReac('S001_S011', surfsys, olhs=[Ca], slhs=[S001], srhs=[S011])
reac11_b = smodel.SReac('S011_S001', surfsys, slhs=[S011], orhs=[Ca], srhs=[S001])

#12 S101 + Ca <=> S111
reac12_f = smodel.SReac('S101_S111', surfsys, olhs=[Ca], slhs=[S101], srhs=[S111])
reac12_b = smodel.SReac('S111_S101', surfsys, slhs=[S111], orhs=[Ca], srhs=[S101])




# Rate constants and their units
reac1_f.kcst = 4e8       # 1/(Ms)
reac1_b.kcst = 52        # 1/s
reac2_f.kcst = 2e5       # 1/(Ms)
reac2_b.kcst = 0.21      # 1/s
reac3_f.kcst = 4e8       # 1/(Ms)
reac3_b.kcst = 377.36    # 1/s
reac4_f.kcst = 2e5       # 1/Ms
reac4_b.kcst = 0.0289    # 1/s
reac5_f.kcst = 4e8       # 1/Ms
reac5_b.kcst = 52        # 1/s
reac6_f.kcst = 2e5       # 1/(Ms)
reac6_b.kcst = 0.21      # 1/s
reac7_f.kcst = 4e8       # 1/(Ms)
reac7_b.kcst = 377.36    # 1/s
reac8_f.kcst = 2e5       # 1/(Ms)
reac8_b.kcst = 0.0289    # 1/s
reac9_f.kcst = 2e7       # 1/(Ms)
reac9_b.kcst = 1.6468    # 1/s
reac10_f.kcst = 2e7      # 1/(Ms)
reac10_b.kcst = 1.6468   # 1/s
reac11_f.kcst = 2e7      # 1/(Ms)
reac11_b.kcst = 1.6468   # 1/s
reac12_f.kcst = 2e7      # 1/(Ms)
reac12_b.kcst = 1.6468   # 1/s

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

