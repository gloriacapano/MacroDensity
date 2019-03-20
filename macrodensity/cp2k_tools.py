from __future__ import print_function
import numpy as np
from itertools import chain
import math


def read_cube_density(FILE, quiet = False):

    print("Reading cube file...")
    with open(FILE, "r") as f:
        _ = f.readline()
        _ = f.readline()
        num_atoms = int(f.readline().split()[0])

        lattice = np.zeros(shape=(3,3))
        tmp = f.readline().split()
        NGX = int(tmp[0])
        lattice[0:1] = float(tmp[1]), float(tmp[2]), float(tmp[3])
        tmp = f.readline().split()
        NGY = int(tmp[0])
        lattice[1:2] = float(tmp[1]), float(tmp[2]), float(tmp[3])
        tmp = f.readline().split()
        NGZ = int(tmp[0])
        lattice[2:3] = float(tmp[1]), float(tmp[2]), float(tmp[3])

        bohr = 0.52917721067
        lattice[0:1]=NGX*lattice[0:1]*bohr
        lattice[1:2]=NGY*lattice[1:2]*bohr
        lattice[2:3]=NGZ*lattice[2:3]*bohr


        atom_type = np.zeros(num_atoms)
        coord = np.zeros(shape=(num_atoms, 3))
        for i in range(num_atoms):
            tmp = f.readline().split()
            atom_type[i] = int(tmp[0])
            coord[i,0] = float(tmp[2])
            coord[i,1] = float(tmp[3])
            coord[i,2] = float(tmp[4])


        print("Reading 3D data...")
        Potential = (f.readline().split()
            for i in range((int(NGZ/6) + (NGZ%6 > 0)) * NGY *NGX))
        Potential = np.fromiter(chain.from_iterable(Potential), float)

    print("Average of the potential = ", np.average(Potential))
    return Potential, NGX, NGY, NGZ, lattice



def density_grid_cube(Density, nx, ny, nz, Volume=1):
    '''Convert the potential
    Args:
       Density: Array of the grid potential of a cube file
       nx,y,z : Number of mesh points in x/y/z
    Returns:
       Potential_grid: the (normalized) quantity on a mesh
    '''
    l = 0
    hartree2eV = 27.211399
    Potential_grid = np.zeros(shape=(nx,ny,nz))
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                Potential_grid[i,j,k] = (Density[l] / Volume) * hartree2eV 
                l = l + 1
    return Potential_grid

