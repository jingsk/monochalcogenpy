from ase.spacegroup.spacegroup import Spacegroup
import numpy as np

class Spacegroup_MX(Spacegroup):
    def __init__(self, sg_no=31, orientation = 'ac'):
        super().__init__(sg_no)
        self.sg_no = sg_no
        if self.sg_no == 31:
            self.ref_orientation = 'cb'
        self.orientation = orientation
        self.sort_index = sort_axis_indices(
            fill_abc(self.orientation),
            fill_abc(self.ref_orientation)
        )
    def get_symop(self):
        rots = self._sorted_rotations()
        trans = self._sorted_translations()
        return [(r,t) for r,t in zip(rots, trans)]
    def _sorted_translations(self):
        _, def_trans = zip(*super().get_symop())
        def_trans = np.array(def_trans)
        trans = def_trans[:,self.sort_index]
        return tuple(trans)
    def _sorted_rotations(self):
        def_rots, _ = zip(*super().get_symop())
        rots = np.array(def_rots).copy()
        #trans = def_trans[:,self.sort_index]
        for i, arr in enumerate(rots):
            diag_elements = np.diag(arr)
            rots[i] = np.diag(diag_elements[self.sort_index])
        return tuple(rots)
        
        
def fill_abc(string):
    assert len(string) ==2
    missing_char = [char for char in 'abc' if char not in string][0]
    return string + missing_char

def sort_axis_indices(orientation, ref_orientation):
    #goal - turn indices currently along abc into current orientation
    #understand this as (the first part gives mapping from 'abc' to ref_orientation.)
    #then the second maps orientation to bc 
    #so in effect we are mapping from orientation to abc obeying order from ref_orientation
    return sort_string_indices('abc',ref_orientation)[sort_string_indices(orientation, 'abc')]

def sort_string_indices(string, ref_string):
    """
    Sorts the indices of string1 based on the sorted order of characters in string2.
    
    Args:
    string: The string whose indices need to be sorted.
    ref_string: The string based on whose sorted order the indices will be sorted.
    
    Returns:
    A list of integers representing the sorted indices of string.
    """
    indices = {c: i for i, c in enumerate(list(string))}
    sort_index = np.array([indices[char] for char in list(ref_string)])
    return sort_index
