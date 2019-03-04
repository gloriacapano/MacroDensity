#! /usr/bin/env python
import macrodensity as md
import math
import numpy as np
import matplotlib.pyplot as plt
from macrodensity.cp2k_tools import read_cube_density


input_file = '../../TEST.cube'
cube_size = [2,2,2]    # This size is in units of mesh points
## origin defines the bottom left point of the cube the "0,0,0" point in fractional coordinates
cube_origin = [0.5,0.5,0.5]
# No need to alter anything after here
#------------------------------------------------------------------
# Get the potential
# This section should not be altered
#------------------------------------------------------------------
pot, NGX, NGY, NGZ, Lattice = read_cube_density(input_file)
print(Lattice)
print(NGX, NGY, NGZ)
vector_a,vector_b,vector_c,av,bv,cv = md.matrix_2_abc(Lattice)
print(vector_a,vector_b,vector_c,av,bv,cv)
resolution_x = vector_a/NGX
resolution_y = vector_b/NGY
resolution_z = vector_c/NGZ
grid_pot = md.density_grid_cube(pot,NGX,NGY,NGZ)
#------------------------------------------------------------------

##------------------------------------------------------------------
# Getting the average potential in a single cube of arbitrary size
##------------------------------------------------------------------
## cube defines the size of the cube in units of mesh points (NGX/Y/Z)
cube = cube_size
## origin defines the bottom left point of the cube the "0,0,0" point in fractional coordinates
origin = cube_origin
## travelled; do not alter this variable
travelled = [0,0,0]
## Uncomment the lines below to do the business
cube_potential, cube_var = md.cube_potential(origin,travelled,cube,grid_pot,NGX,NGY,NGZ)
print "Potential            Variance"
print "--------------------------------"
print cube_potential,"   ", cube_var
##------------------------------------------------------------------
##------------------------------------------------------------------
## Uncomment the lines below to do the business
dim = [1,2,3,4,5,6,7,8,9,10,20,40,60,80,100]
print ("Dimension   Potential   Variance")
print ("--------------------------------")
for d in dim:
    cube = [d,d,d]
    cube_potential, cube_var = md.cube_potential(cube_origin,travelled,cube,grid_pot,NGX,NGY,NGZ)
    print(" %3i     %10.4f   %10.6f"%(d,cube_potential,cube_var))
##------------------------------------------------------------------
##------------------------------------------------------------------

