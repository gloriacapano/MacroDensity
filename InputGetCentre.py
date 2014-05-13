import PotentialModule as pot
import math
import numpy as np
import matplotlib.pyplot as plt





#------------------------------------------------------------------
#   READING
# Get the two potentials and change them to a planar average.
# This section should not be altered
#------------------------------------------------------------------
## SLAB
vasp_pot, NGX, NGY, NGZ, Lattice = pot.read_vasp_density('CHGCAR.slab')
mag_a,mag_b,mag_c,vec_a,vec_b,vec_c = pot.matrix_2_abc(Lattice)
resolution_x = mag_a/NGX
resolution_y = mag_b/NGY
resolution_z = mag_c/NGZ
Volume = pot.get_volume(vec_a,vec_b,vec_c)
grid_pot_slab, electrons_slab = pot.density_2_grid(vasp_pot,NGX,NGY,NGZ,True,Volume)
# Save the lattce vectors for use later
Vector_A = [vec_a,vec_b,vec_c]
##----------------------------------------------------------------------------------
## CONVERT TO PLANAR DENSITIES
##----------------------------------------------------------------------------------
slab = pot.planar_average(grid_pot_slab,NGX,NGY,NGZ)
#----------------------------------------------------------------------------------
## BULK
vasp_pot, NGX, NGY, NGZ, Lattice = pot.read_vasp_density('CHGCAR.bulk')
mag_a,mag_b,mag_c,vec_a,vec_b,vec_c = pot.matrix_2_abc(Lattice)
resolution_x = mag_a/NGX
resolution_y = mag_b/NGY
resolution_z = mag_c/NGZ
Volume = pot.get_volume(vec_a,vec_b,vec_c)
grid_pot_bulk, electrons_bulk = pot.density_2_grid(vasp_pot,NGX,NGY,NGZ,True,Volume)
# Save the lattce vectors for use later
Vector_B = [vec_a,vec_b,vec_c]
##----------------------------------------------------------------------------------
## CONVERT TO PLANAR DENSITIES
##----------------------------------------------------------------------------------
bulk = pot.planar_average(grid_pot_bulk,NGX,NGY,NGZ)
#----------------------------------------------------------------------------------
# FINISHED READING
#----------------------------------------------------------------------------------
# MOVE THE SLAB TO AVOID PBC
#----------------------------------------------------------------------------------
slab, bulk = pot.matched_spline_generate(slab, bulk, Vector_A[2], Vector_B[0])
slab = pot.translate_grid(slab,10.0,True,Vector_A[2],0.0)
#----------------------------------------------------------------------------------
# TEST PLOT TO ENSURE THAT THE DENSITIES ARE CLOSE ENOUGH
#----------------------------------------------------------------------------------
#plt.plot(bulk[:,0],bulk[:,1],'o')
#plt.plot(slab[:,0],slab[:,1],'*')
#plt.show()

#----------------------------------------------------------------------------------
# GET RATIOS
#----------------------------------------------------------------------------------
elect_ratio = (electrons_slab/electrons_bulk)
print "RATIO: ",elect_ratio
extended = pot.symmetric_extend_potential(bulk,round(elect_ratio,8)) # Rounding the ratio gives a much smoother result.
#----------------------------------------------------------------------------------
# TEMP SPLINE TO VERY HIGH RESOLUTION
# This is to ensure accurate matching of centres
#----------------------------------------------------------------------------------
print "Generating high-res splines"
extended = pot.c_spline_generate(extended,0.005)
slab = pot.c_spline_generate(slab,0.005)
#----------------------------------------------------------------------------------
# OVERLAY THE TWO DENSITIES
#-----------------------------------------------------------------------------------
# A) Centre of the extended bulk
if len(extended) % 2 != 0:
    centre_of_bulk = extended[len(extended)/2,0]
else:
    centre_of_bulk = (extended[len(extended)/2,0]+extended[len(extended)/2+1,0])/2 
print "CENTRE OF BULK: ", centre_of_bulk
# B) Centre of the slab
centre_of_charge = pot.centre_of_charge(slab,tol=0.000005)
# The index of that point:
centre_point = int(len(slab)/((slab[len(slab)-1,0]-slab[0,0])/centre_of_charge))
# C) Matching
bulk = pot.translate_grid(extended, centre_of_charge-centre_of_bulk,False,Vector_B[2])
centre_of_bulk = extended[len(extended)/2,0]
print "TRANSLATED CENTRE OF BULK: ", centre_of_bulk
# D) Refined matching of the centres, based on matching extrema
new_offset = pot.match_extrema(bulk,slab,len(bulk)/2,centre_point)
# Spline to a lower resolution for the rest of the analysis
print "Reducing resolutions"
print "Original bulk length: ",len(bulk)
bulk = pot.spline_generate(bulk,10)
slab = pot.spline_generate(slab,10)
print "Modified bulk length: ",len(bulk)
#slab = pot.translate_grid(slab,0.034,False,Vector_A[2],0.0)
bulk = pot.translate_grid(bulk,new_offset,True,elect_ratio*Vector_B[0],0.0)
plt.plot(bulk[:,0],bulk[:,1],'o')
plt.plot(slab[:,0],slab[:,1],'*')
plt.show()

#-----------------------------------------------------------------------------------
# ADD THE VACUUM SECTION OUTSIDE THE BULK
#-----------------------------------------------------------------------------------
bulk = pot.bulk_vac(bulk, slab)
#-----------------------------------------------------------------------------------
#SPLINE TO THE SAME RESOLUTION
#-----------------------------------------------------------------------------------
Res_ratio = len(bulk)/len(slab)
bulk = pot.spline_generate(bulk,Res_ratio)
#-----------------------------------------------------------------------------------
# GET THE DIFFERENCE
#-----------------------------------------------------------------------------------

print "Upping the resolution for difference"
print "Original bulk length: ",len(bulk)
bulk = pot.spline_generate(bulk,0.1)
slab = pot.spline_generate(slab,0.1)
print "Modified bulk length: ",len(bulk)
#plt.plot(slab[:,0],slab[:,1],'o')
#plt.plot(bulk[:,0],bulk[:,1],'*')
#plt.plot(bulk[:,0],bulk[:,1],'*')
plt.show()
difference = pot.subs_potentials(slab,bulk,tol=1.0e-8)
difference = pot.spline_generate(difference,100)
plt.plot(difference[:,0],difference[:,1])
plt.show()
np.savetxt('ChargeDifference.dat',difference)
