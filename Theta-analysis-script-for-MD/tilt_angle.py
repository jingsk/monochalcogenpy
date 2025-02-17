import numpy as np
from itertools import permutations

supported_projection = [v+w for v,w in permutations('abc', 2)]

def _projected_vector_angle(v, proj ='ab'):
    proj = proj.lower()
    ordered_idx = list(proj) +[axis for axis in 'abc' if axis not in proj]
    opp_idx, adj_idx, flatted_idx = ['abc'.find(i) for i in ordered_idx]
    v[flatted_idx] = 0
    #cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arctan(v[opp_idx]/v[adj_idx])
    return np.degrees(angle)

def tilt_angle(atoms, proj):
    '''
    proj is the plane to project v,w and the string 
    has the following order: opposite, adjacent
    '''
    if proj not in supported_projection:
        raise ValueError(f"Projection {proj} not supported.")
    #return projected_vector_angle(*[atoms.get_positions()[i] for i in [2,5]])