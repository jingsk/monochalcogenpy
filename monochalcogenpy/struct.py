#so empty
from ase import Atoms
import numpy as np
import warnings

class SpacegroupError(Exception):
    """Base exception for the spacegroup module."""

class SpacegroupValueError(SpacegroupError):
    """Raised when arguments have invalid value."""


def equivalent_sites(basis, spg_dataset,symprec = 1e-4):
    """Create an Atoms instance for a conventional unit cell of a
    space group. Originally from ase modified by Jing T.

    Parameters:

    basis : list of scaled coordinates
        Positions of the unique sites corresponding to symbols given
        either as scaled positions or through an atoms instance.  Not
        needed if *symbols* is a sequence of Atom objects or an Atoms
        object.
        symbols : str | sequence of str | sequence of Atom | Atoms
        Element symbols of the unique sites.  Can either be a string
        formula or a sequence of element symbols. E.g. ('Na', 'Cl')
        and 'NaCl' are equivalent.  Can also be given as a sequence of
        Atom objects or an Atoms object.
    spg_dataset : spglib dataset obtained from ase.spacegroup.symmetrize.check_symmetry
    """
    #from ase.spacegroup.symmetrize import check_symmetry
    #s_ds = check_symmetry(atoms)
    rotations = spg_dataset.rotations
    translations = spg_dataset.translations
    
    def find_orbit(point: np.ndarray) -> np.ndarray:
        """Find crystallographic orbit of the given point."""
        candidates = ((rotations @ point) + translations % 1.0) % 1.0
        orbit = [candidates[0]]
        for member in candidates[1:]:
            diff = member - orbit
            diff -= np.rint(diff)
            if not np.any(np.all(np.abs(diff) < symprec, axis=1)):
                orbit.append(member)
        return np.array(orbit)
    
    onduplicates = 'keep'
    orbits = []
    for kind, pos in enumerate(basis):
        for i, (kind0, positions0) in enumerate(orbits):
            diff = pos - positions0
            diff -= np.rint(diff)
            if np.any(np.all(np.abs(diff) < symprec, axis=1)):
                if onduplicates == 'keep':
                    pass
                elif onduplicates == 'replace':
                    orbits[i] = (kind, positions0)
                elif onduplicates == 'warn':
                    warnings.warn(
                        'scaled_positions %d and %d are equivalent' %
                        (kind0, kind))
                elif onduplicates == 'error':
                    raise SpacegroupValueError(
                        'scaled_positions %d and %d are equivalent' %
                        (kind0, kind))
                break
        else:
            orbits.append((kind, find_orbit(pos)))
    
    kinds = []
    sites = []
    for kind, orbit in orbits:
        kinds.extend(len(orbit) * [kind])
        sites.append(orbit)
    
    return np.concatenate(sites, axis=0), kinds



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