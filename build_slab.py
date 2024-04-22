from ase.io import lammpsdata as lmp
from ase.build import fcc111

slab = fcc111('Pt', size=(12,14,6), vacuum=15.0, orthogonal= True)

lmp.write_lammps_data('pt-standard-slab.lmp', atoms = slab, units = 'metal', atom_style = 'charge')