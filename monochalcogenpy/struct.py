#so empty
from ase import Atoms
import numpy as np

def unit_cell(a, b, c, use_symm=False):
    cell = np.diag([a,b,c])
    #define Se-Se height 
    h = np.sqrt(3 - (c/a)**2)/2 * a /b
    #angle relative of Se-Ge vector to vdw direction along the ac plane
    theta = np.arctan(np.sqrt(2)) - np.arctan(2*h*b/c)
    #displacement 
    z_Ge = 1/np.sqrt(2) * a * np.cos(theta) / b
    y_Ge = 1/np.sqrt(2) * a * np.sin(theta) / c
    if not use_symm:
        pos = np.array([
            [0, 0.5 - h / 2+z_Ge, 0 + y_Ge], #Ge1
            [0.5, 0.5 + h / 2-z_Ge, 0.5 + y_Ge], #Ge2
            [0, 0.5 - h / 2, 0],  #Se1
            [0.5, 0.5 + h / 2, 0.5],  #Se2
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
            [0, 0.5 - h / 2+z_Ge, 0 + y_Ge], #Ge
            [0, 0.5 - h / 2, 0]              #Se
        ]) 
        atoms = crystal(
        ('Ge', 'Se'), 
        basis=basis, 
        spacegroup = 31, 
        cell=cell)
    return atoms