#so empty
from ase import Atoms
import numpy as np
from monochalcogenpy.crystal import crystal
from monochalcogenpy.utils import fill_abc, sort_axis_indices
from monochalcogenpy.spacegroup import Spacegroup_MX

def unit_cell(a, b, c, orientation ='ac', use_symm=False):
    """
    Parameters
    ----------
    a: Float
        vector length along a
    b: Float
        vector length along b
    c: Float
        vector length along c 
    orientation: Str
        orientation of the armchair direction followed by 
        the direction normal to the layer. For example 
        orientation ='ac' means armchair is along a and vdw along c. 
        orientation ='bc' means armchair is along b and vdw along c.
    use_symm: Bool
        whether to use symmetry-based 2-atom basis definition or to 
        explicitly define 4-atom scaled position. Either gives the same
        answer but one retains spacegroup symmetry.

    Scaled atomic positions are used to create a unit cell of monolayer monochalcogen.
    The implementaiton assumes all spacegroup 31 symmetry operations. Relative positions
    are obtained from rigid polyhedral approximation. Two implementations are coded and 
    they give the same answer but use_symm returns Atoms with spacegroup symmetry. 
    
    Currently ref_orientation is written to be a1= armchair direction and a2=zigzag but 
    output atomic positions are rotated to specified orientation.
    """
    ref_orientation = 'ac'
    #currently only ac orientation is supported
    cell = np.diag([a,b,c])
    #define a1 > a2
    [a2, a1] = np.sort([a,b])
    #define Se-Se height 
    h = np.sqrt(3 - (a1/a2)**2)/2 * a2 /c
    #angle relative of Se-Ge vector to vdw direction along the ac plane
    theta = np.arctan(np.sqrt(2)) - np.arctan(2*h*c/a1)
    #displacement 
    z_Ge = 2.56 * np.cos(theta) / c
    x_Ge = 2.56 * np.sin(theta) / a1
    if not use_symm:
        pos = np.array([
            [0   + x_Ge, 0,   0.5 - h / 2+z_Ge, ], #Ge1
            [0.5 + x_Ge, 0.5, 0.5 + h / 2-z_Ge], #Ge2
            [0,          0,   0.5 - h / 2],  #Se1
            [0.5,        0.5, 0.5 + h / 2],  #Se2
            ])
        pos = pos_by_orientation(pos, orientation, ref_orientation)
        atoms = Atoms(
            'Ge2Se2',
            scaled_positions=pos,
            cell = cell,
            pbc=[True,True,True]
            )
    else:
        basis = np.array([
            [0 + x_Ge, 0, 0.5 - h / 2+z_Ge], #Ge
            [0,        0, 0.5 - h / 2]       #Se
        ])
        basis = pos = pos_by_orientation(basis, orientation, ref_orientation)
        sg = Spacegroup_MX(sg_no=31, orientation=orientation)
        atoms = crystal(
        ('Ge', 'Se'), 
        basis=basis, 
        spacegroup = sg, 
        cell=cell)
    return atoms


def unit_cell_bulk(a, b, c, registry=[0., 0.], orientation ='ac', use_symm=True):
    """
    Parameters
    ----------
    a: Float
        vector length along a
    b: Float
        vector length along b
    c: Float
        vector length along c 
    registry: List
        Scaled position shifts (interlayer registry) along the armchair 
        and zigzag direction respectively.
    orientation: Str
        orientation of the armchair direction followed by 
        the direction normal to the layer. For example 
        orientation ='ac' means armchair is along a and vdw along c. 
        orientation ='bc' means armchair is along b and vdw along c.
    use_symm: Bool
        whether to use symmetry-based 2-atom basis definition or to 
        explicitly define 8-atom scaled position. Either gives the same
        answer but one retains spacegroup symmetry.

    Scaled atomic positions are used to create a unit cell of bulk monochalcogen.
    The implementaiton assumes all spacegroup 62 symmetry operations. Relative positions
    are obtained from rigid polyhedral approximation. Interlayer registry are implemented as half the shift. 
    Two implementations are coded and they give the same answer but use_symm returns Atoms with spacegroup symmetry. 
    
    Currently ref_orientation is written to be a1= armchair direction and a2=zigzag but 
    output atomic positions are rotated to specified orientation.
    """
    ref_orientation = 'ac'
    #currently only ac orientation is supported
    cell = np.diag([a,b,c])
    #define a1 > a2
    [a2, a1] = np.sort([a,b])
    #define Se-Se height 
    h = np.sqrt(3 - (a1/a2)**2)/2 * a2 /c
    #angle relative of Se-Ge vector to vdw direction along the ac plane
    theta = np.arctan(np.sqrt(2)) - np.arctan(2*h*c/a1)
    #displacement 
    z_Ge = 2.56 * np.cos(theta) / c
    x_Ge = 2.56 * np.sin(theta) / a1

    if not use_symm:
        pos = np.array([
            [0.0 + x_Ge, 0.25, 0.25 - h / 2+z_Ge], #Ge1
            [0.5 + x_Ge, 0.75, 0.25 + h / 2-z_Ge], #Ge2
            [0.0 - x_Ge, 0.75, 0.75 - h / 2+z_Ge], #Ge3
            [0.5 - x_Ge, 0.25, 0.75 + h / 2-z_Ge], #Ge4
            [0.5,        0.75, 0.25 - h / 2],      #Se1
            [0.0,        0.25, 0.25 + h / 2],      #Se2
            [0.5,        0.25, 0.75 + h / 2],      #Se3
            [0.0,        0.75, 0.75 - h / 2],      #Se4            
            ])
        pos = pos_by_orientation(pos, orientation, ref_orientation)
        atoms = Atoms(
            'Ge4Se4',
            scaled_positions=pos,
            cell = cell,
            pbc=[True,True,True]
            )
    else:
        basis = np.array([
            [0 + x_Ge, 0.25, 0.25 + h / 2-z_Ge], #Ge
            [0.5,      0.75, 0.25 - h / 2]       #Se
        ])
        basis += [registry[0]/2, registry[1]/2, 0]
        basis = pos_by_orientation(basis, orientation, ref_orientation)
        sg = Spacegroup_MX(sg_no=62, orientation=orientation)
        atoms = crystal(
        ('Ge', 'Se'), 
        basis=basis, 
        spacegroup = sg, 
        cell=cell)
    return atoms

def pos_by_orientation(pos, orientation, ref_orientation):
    """
    Rotate positions currently specified in assumed ref_orientation 
    to specified (desired) orientation.
    """
    sort_idx = sort_axis_indices(fill_abc(orientation), fill_abc(ref_orientation))
    return pos[:,sort_idx]