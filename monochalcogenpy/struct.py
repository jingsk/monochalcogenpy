#so empty
from ase import Atoms
import numpy as np

def unit_cell(a, b, c, use_symm=False):
    #replacing #b with c
    #replacing #c with a
    #replacing #a with b
    cell = np.diag([a,b,c])
    #define Se-Se height 
    h = np.sqrt(3 - (a/b)**2)/2 * b /c
    #angle relative of Se-Ge vector to vdw direction along the ac plane
    theta = np.arctan(np.sqrt(2)) - np.arctan(2*h*c/a)
    #displacement 
    z_Ge = 1/np.sqrt(2) * b * np.cos(theta) / c
    x_Ge = 1/np.sqrt(2) * b * np.sin(theta) / a
    if not use_symm:
        pos = np.array([
            [0   + x_Ge, 0,   0.5 - h / 2+z_Ge, ], #Ge1
            [0.5 + x_Ge, 0.5, 0.5 + h / 2-z_Ge], #Ge2
            [0,          0,   0.5 - h / 2],  #Se1
            [0.5,        0.5, 0.5 + h / 2],  #Se2
            ])
        atoms = Atoms(
            'Ge2Se2',
            scaled_positions=pos,
            cell = cell,
            pbc=[True,True,True]
            )
    else:
        from ase.spacegroup import crystal
        basis = np.array([
            [0 + x_Ge, 0, 0.5 - h / 2+z_Ge], #Ge
            [0,        0, 0.5 - h / 2]              #Se
        ]) 
        atoms = crystal(
        ('Ge', 'Se'), 
        basis=basis, 
        spacegroup = 31, 
        cell=cell)
    return atoms