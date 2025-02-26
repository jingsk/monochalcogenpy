#so empty
from ase import Atoms
import numpy as np

def unit_cell(a, b, c, method):
    cell = np.diag([a,b,c])
    #define Se-Se height 
    h = np.sqrt(3 - (c/a)**2)/2 * a /b
    #angle relative of Se-Ge vector to vdw direction along the ac plane
    theta = np.arctan(np.sqrt(2)) - np.arctan(2*h*b/c)
    #displacement 
    z_Ge = 1/np.sqrt(2) * a * np.cos(theta) / b
    y_Ge = 1/np.sqrt(2) * a * np.sin(theta) / c
    if method == 'pos':
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
        return atoms