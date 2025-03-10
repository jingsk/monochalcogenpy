from ase.spacegroup.spacegroup import Spacegroup
import numpy as np
from monochalcogenpy.utils import fill_abc, sort_axis_indices

class Spacegroup_MX(Spacegroup):
    def __init__(self, sg_no=31, orientation = 'ac'):
        super().__init__(sg_no)
        self.sg_no = sg_no
        if self.sg_no == 31:
            self.ref_orientation = 'cb'
        if self.sg_no == 62:
            self.ref_orientation = 'ca'
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