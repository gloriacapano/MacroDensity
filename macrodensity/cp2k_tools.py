from __future__ import print_function
import numpy as np
from itertools import chain



def read_cube_density(FILE, use_pandas = None, quiet = False):
    if use_pandas:
        from pandas import read_table as pandas_read_table
    elif use_pandas is None:
        try:
            from pandas import read_table as pandas_read_table
            use_pandas = True
        except ImportError:
            use_pandas = False

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


        atom_type = np.zeros(num_atoms)
        coord = np.zeros(shape=(num_atoms, 3))
        for i in range(num_atoms):
            tmp = f.readline().split()
            atom_type = int(tmp[0])
            coord[i,0] = float(tmp[1])
            coord[i,1] = float(tmp[2])
            coord[i,2] = float(tmp[3])

        if use_pandas:
            print("Reading 3D data using Pandas...")
            skiprows = 10 + num_atoms
            readrows = int(math.ceil(NGX * NGY * NGZ / 6))

            dat = pandas_read_table(FILE, delim_whitespace=True,
                                    skiprows=skiprows, header=None,
                                    nrows=readrows)
            Potential = dat.iloc[:readrows, :6].values.flatten()
            remainder = (NGX * NGY * NGZ) % 6
            if remainder > 0:
                Potential = Potential[:(-6 + remainder)]

        else:
            print("Reading 3D data...")
            Potential = (f.readline().split()
                             for i in range(int(math.ceil(NGX * NGY * NGZ / 6))))
            Potential = np.fromiter(chain.from_iterable(Potential), float)

    _print_boom(quiet=quiet)
    if not quiet:
        print("Average of the potential = ", numpy.average(Potential))

    return Potential, NGX, NGY, NGZ, lattice

