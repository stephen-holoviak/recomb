from ase import Atom, Atoms
import ase.io.lammpsdata as lammps

# read the initial structure
slab = lammps.read_lammps_data(file='pt-asym-init_water.lammps', units='metal', style='charge')

# make delete the top two Pt layers from one slab to create an asymmetric simulation
del slab[[atom.index for atom in slab if atom.number == 3 and atom.z > 119]]

# ---------- find the z level of each layer in the electrodes
metal_type = 3
unique_Zs = []

for atom in slab:
    if atom.number == metal_type and round(atom.z,0) not in unique_Zs:
        unique_Zs.append(round(atom.z,0))
# print(unique_Zs)
unique_Zs.sort()
# ---------- Change the id number of each layer to proper group
for atom in slab:
    if atom.number == metal_type:
        if round(atom.z,0) == unique_Zs[0]: #anode fixed
            atom.number = 3
        if round(atom.z,0) in unique_Zs[1:3]: #anode thermo
            atom.number = 4
        if round(atom.z,0) in unique_Zs[3:8]: #anode work
            atom.number = 5
        if round(atom.z,0) in unique_Zs[8:11]: #cathode work
            atom.number = 6
        if round(atom.z,0) in unique_Zs[11:13]: #cathode thermo
            atom.number = 7
        if round(atom.z,0) == unique_Zs[13:]: #cathode fix
            atom.number = 8


lammps.write_lammps_data('asym_water.input', atoms = slab, units = 'metal', atom_style = 'charge')

