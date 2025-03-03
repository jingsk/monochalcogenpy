#so empty
from ase import Atoms
import numpy as np
from monochalcogenpy.crystal import crystal
from monochalcogenpy.spacegroup import Spacegroup_MX

def unit_cell(a, b, c, orientation ='ac', use_symm=False):
    #currently only ac orientation is supported
    cell = np.diag([a,b,c])
    #define Se-Se height 
    h = np.sqrt(3 - (a/b)**2)/2 * b /c
    #angle relative of Se-Ge vector to vdw direction along the ac plane
    theta = np.arctan(np.sqrt(2)) - np.arctan(2*h*c/a)
    #displacement 
    z_Ge = 2.56 * np.cos(theta) / c
    x_Ge = 2.56 * np.sin(theta) / a
    if not use_symm:
        if orientation == 'ac':
            pos = np.array([
                [0   + x_Ge, 0,   0.5 - h / 2+z_Ge, ], #Ge1
                [0.5 + x_Ge, 0.5, 0.5 + h / 2-z_Ge], #Ge2
                [0,          0,   0.5 - h / 2],  #Se1
                [0.5,        0.5, 0.5 + h / 2],  #Se2
                ])
        if orientation == 'bc':
            pos = np.array([
                [0,   0 + x_Ge,   0.5 - h / 2+z_Ge, ], #Ge1
                [0.5, 0.5 + x_Ge, 0.5 + h / 2-z_Ge], #Ge2
                [0,   0,          0.5 - h / 2],  #Se1
                [0.5, 0.5,        0.5 + h / 2],  #Se2
                ])
        atoms = Atoms(
            'Ge2Se2',
            scaled_positions=pos,
            cell = cell,
            pbc=[True,True,True]
            )
    else:
        if orientation == 'ac':
            basis = np.array([
                [0 + x_Ge, 0, 0.5 - h / 2+z_Ge], #Ge
                [0,        0, 0.5 - h / 2]       #Se
            ])
        if orientation == 'bc':
            basis = np.array([
                [0, 0  + x_Ge, 0.5 - h / 2+z_Ge], #Ge
                [0, 0        , 0.5 - h / 2]       #Se
            ])
        sg = Spacegroup_MX(sg_no=31, orientation=orientation)
        atoms = crystal(
        ('Ge', 'Se'), 
        basis=basis, 
        spacegroup = sg, 
        cell=cell)
    return atoms