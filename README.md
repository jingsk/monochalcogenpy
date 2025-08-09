## Monochalcogenpy
Chemical structure generation and analysis package following proper spacegroup symmetry. 

## Geometric description of monochalcogenide
This contruction of monochalcogenide structure describes each layer as an intermediate structure between h-BN, a flat graphene-like layer, and black phosphorus. As the ratio between the two in-plane lattice parameters decrease, Se atoms are pushed out of plane creating truss-like scaffold with Ge atoms as the tip of trigonal pyramids pointing up and down alternately. 

*images coming*



## How to build structure

To build a bulk GeSe with lattice vectors (a,b,c) with the armchair direction along $a$ and out of plane along direction $c$ with spacegroup 31 symmetry enforcement.
` atoms = unit_cell(a, b, c, orientation ='ac', use_symm=True)`
To build a monolayer GeSe with lattice vectors (a,b,c) with the armchair direction along $c$ and out of plane direction along $a$ with spacegroup 62 symmetry enforcement.
` atoms = unit_cell_bulk(a, b, c, orientation ='ca', use_symm=True)`

Measure MXM tilt angle with the armchair direction along $a$ and out of plane along direction $c$ as defined in 10.1103/PhysRevB.97.024110:
`theta_deg = tilt_angle(atoms, proj = 'ac', thres=0.9, absolute=True)`
