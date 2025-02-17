import numpy as np
from itertools import permutations
from matscipy.neighbours import neighbour_list

supported_projection = [v+w for v,w in permutations('abc', 2)]

def _projected_vector_angle(v, proj ='ab'):
    proj = proj.lower()
    ordered_idx = list(proj) +[axis for axis in 'abc' if axis not in proj]
    opp_idx, adj_idx, flatted_idx = ['abc'.find(i) for i in ordered_idx]
    v[flatted_idx] = 0
    #cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.abs(np.arctan(v[opp_idx]/v[adj_idx]))
    return np.degrees(angle)

def vec_dir_check(v, dir ='c', thres = 0.8):
    dir = dir.lower()
    dir_idx = 'abc'.find(dir)
    ref_vec = [0,0,0]
    ref_vec[dir_idx] = 1
    return v@ref_vec/np.linalg.norm(v) > thres

def tilt_angle(atoms, proj):
    '''
    proj is the plane to project v,w and the string 
    has the following order: opposite, adjacent
    '''
    if proj not in supported_projection:
        raise ValueError(f"Projection {proj} not supported.")
    src, dst, dis_vec = neighbour_list(
        'ijD',
        atoms=atoms,
        cutoff={('Ge', 'Se'): 3}
    )
    Ge_idx = [atom.index for atom in atoms if atom.symbol=='Ge']
    tilt_angle = {}
    for idx in Ge_idx:
        for D in dis_vec[src==idx]:
            if vec_dir_check(D, dir =proj[1], thres = 0.9):
                try: 
                    tilt_angle[idx]
                    print(f'Duplicate found at idx = {idx}.')
                    print(f'dis_vec = {dis_vec[src == idx]}.')
                    tilt_angle[idx] = tilt_angle[idx]+[_projected_vector_angle(D, proj)]
                except KeyError:
                    tilt_angle[idx] = [_projected_vector_angle(D, proj)]
    return tilt_angle
                