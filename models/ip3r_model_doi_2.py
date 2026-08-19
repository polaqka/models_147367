# ip3r_model.py
#
# Katri Hituri
# IP3R model by Doi et al. 2005
###########

# Import packages

import steps.model as smodel
import steps.geom as sgeom



# *MODEL*
#

mdl = smodel.Model()
volsys = smodel.Volsys('vsys', mdl) # Volume system
surfsys = smodel.Surfsys('ssys', mdl) # Surface system


# CHEMICAL SPECIES 
Ca = smodel.Spec('Ca', mdl)				# Calcium
IP3 = smodel.Spec('IP3', mdl)				# IP3

# IP3 receptor states
R = smodel.Spec('R', mdl)					# IP3 receptor with no bound ligands
RIP3 = smodel.Spec('RIP3', mdl)				# bound IP3 
Ropen = smodel.Spec('Ropen', mdl)			# bound IP3 and Ca (open)
RCa = smodel.Spec('RCa', mdl)				# 1 bound Ca to inactivation site
RCa2 = smodel.Spec('RCa2', mdl)				# 2 bound Ca to inactivation sites
RCa3 = smodel.Spec('RCa3', mdl)				# 3 bound Ca to inactivation sites
RCa4 = smodel.Spec('RCa4', mdl)				# 4 bound Ca to inactivation sites

# NEW
RCa_IP3 = smodel.Spec('RCa_IP3', mdl) # bound IP3 and 1 Ca to inactivation site (closed)


# REACTIONS


# RIP3(s) + Ca <=> Ropen
RIP3_bind_Ca_f = smodel.SReac('RIP3_bind_Ca_f', surfsys, \
			      olhs=[Ca], slhs=[RIP3], srhs = [Ropen])
RIP3_bind_Ca_b = smodel.SReac('RIP3_bind_Ca_b', surfsys, \
			      slhs=[Ropen], orhs=[Ca], srhs=[RIP3])

# R(s) + IP3 <=> RIP3(s)
R_bind_IP3_f = smodel.SReac('R_bind_IP3_f', surfsys, \
			    olhs=[IP3], slhs=[R], srhs=[RIP3])
R_bind_IP3_b = smodel.SReac('R_bind_IP3_b', surfsys, \
			    slhs=[RIP3], orhs=[IP3], srhs=[R])

# R(s) + Ca <=> RCa(s)
R_bind_Ca_f = smodel.SReac('R_bind_Ca_f', surfsys, \
			   olhs=[Ca], slhs=[R], srhs=[RCa])
R_bind_Ca_b = smodel.SReac('R_bind_Ca_b', surfsys, \
			   slhs=[RCa], orhs=[Ca], srhs=[R])

# RCa(s) + Ca <=> RCa2(s)
RCa_bind_Ca_f = smodel.SReac('RCa_bind_Ca_f', surfsys, \
			     olhs=[Ca], slhs=[RCa], srhs = [RCa2])
RCa_bind_Ca_b = smodel.SReac('RCa_bind_Ca_b', surfsys, \
			     slhs=[RCa2], orhs=[Ca], srhs=[RCa])

# RCa2(s) + Ca <=> RCa3(s)
RCa2_bind_Ca_f = smodel.SReac('RCa2_bind_Ca_f', surfsys, \
			      olhs=[Ca], slhs= [RCa2], srhs = [RCa3])
RCa2_bind_Ca_b = smodel.SReac('RCa2_bind_Ca_b', surfsys, \
			      slhs=[RCa3], orhs=[Ca], srhs= [RCa2])

# RCa3(s) + Ca <=> RCa4(s)
RCa3_bind_Ca_f = smodel.SReac('RCa3_bind_ca_f', surfsys, \
			      olhs=[Ca], slhs=[RCa3], srhs=[RCa4])
RCa3_bind_Ca_b = smodel.SReac('RCa3_bind_ca_b', surfsys, \
			      slhs=[RCa4], orhs=[Ca], srhs=[RCa3])

# NEW
'''
# RCa(s) + IP3 <=> Ropen
RCa_bind_IP3_f = smodel.SReac('RCa_bind_IP3_f', surfsys, \
			    olhs=[IP3], slhs=[RCa], srhs=[Ropen])
RCa_bind_IP3_b = smodel.SReac('RCa_bind_IP3_b', surfsys, \
			    slhs=[Ropen], orhs=[IP3], srhs=[RCa])
'''

# RCa(s) + IP3 <=> RCa_IP3(s)
RCa_bind_IP3_f = smodel.SReac('RCa_bind_IP3_f', surfsys, \
olhs=[IP3], slhs=[RCa], srhs=[RCa_IP3])
RCa_bind_IP3_b = smodel.SReac('RCa_bind_IP3_b', surfsys, \
slhs=[RCa_IP3], orhs=[IP3], srhs=[RCa])

# RIP3(s) + Ca <=> RCa_IP3(s)
RIP3_bind_Ca_inact_f = smodel.SReac('RIP3_bind_Ca_inact_f', surfsys, \
olhs=[Ca], slhs=[RIP3], srhs=[RCa_IP3])
RIP3_bind_Ca_inact_b = smodel.SReac('RIP3_bind_Ca_inact_b', surfsys, \
slhs=[RCa_IP3], orhs=[Ca], srhs=[RIP3])


# the reaction constants
R_bind_IP3_f.kcst =  1000e6 # 1/(Ms)
R_bind_IP3_b.kcst = 25800 # 1/s
RIP3_bind_Ca_f.kcst = 8000e6	# 1/(Ms)
RIP3_bind_Ca_b.kcst = 2000	# 1/s
R_bind_Ca_f.kcst = 8.889e6  # 1/(Ms)
R_bind_Ca_b.kcst = 5		# 1/s
RCa_bind_Ca_f.kcst = 20e6 	# 1/(Ms)
RCa_bind_Ca_b.kcst = 10 	# 1/s
RCa2_bind_Ca_f.kcst = 40e6 	# 1/(Ms)
RCa2_bind_Ca_b.kcst = 15  	# 1/s
RCa3_bind_Ca_f.kcst = 60e6  	# 1/(Ms)
RCa3_bind_Ca_b.kcst = 20  	# 1/s

# NEW
'''
RCa_bind_IP3_f.kcst =  1000e6  #1000e6 # 1/(Ms)
RCa_bind_IP3_b.kcst = 11466.81 #11466.81 # 1/s
'''

# Kinetyka wiązania IP3 taka sama jak dla pustego receptora (R_bind_IP3)
RCa_bind_IP3_f.kcst = 1000e6 # 1/(Ms)
RCa_bind_IP3_b.kcst = 25800 # 1/s

# Kinetyka wiązania Ca inaktywującego taka sama jak dla pustego receptora (R_bind_Ca)
RIP3_bind_Ca_inact_f.kcst = 8.889e6 # 1/(Ms)
RIP3_bind_Ca_inact_b.kcst = 5 # 1/s
          
    

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

