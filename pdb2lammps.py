from ase.io.proteindatabank import read_proteindatabank
from ase.io.lammpsdata import write_lammps_data
import numpy as np

in_atoms = read_proteindatabank(file='electrolyte.pdb')
in_atoms.set_cell([[33.0, 0.0, 0.0],[0.0, 33.0, 0.0],[0.0, 0.0, 72.0]])
write_lammps_data(file='defect-cu_water.input', atoms=in_atoms, velocities=False, units='metal', atom_style='charge')
