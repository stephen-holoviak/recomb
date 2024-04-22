# Uses a grid to find neighboring atoms for a set of test atoms, should be faster than the old version
# Stephen Holoviak 04/2023
import read_TS as ts
from plot_rdf import pbc_distance
from tqdm import tqdm
import pandas as pd
import os



def group_atoms(data, cell, r=5,  pickle=''):
    groups = {}
    x_max, y_max, z_max = int(cell[0]//r), int(cell[1]//r), int(cell[2]//r)

    for x in range(0, x_max+1):
        for y in range(0, y_max+1):
            for z in range(0, z_max+1):
                cube = (x, y, z)
                if cube not in groups:
                    groups[cube]=[]
    # Loop through each point and assign it to a group based on its position in space
   
    for index, row in data.iterrows():

        # Calculate the indices of the cubic region that the point belongs to
        x_index = int(row['X'] // r)
        if x_index < 0:
            x_index = 0
        y_index = int(row['Y'] // r)
        if y_index < 0:
            y_index = 0
        z_index = int(row['Z'] // r)
        if z_index < 0:
            z_index = 0

        # Create a tuple representing the cubic region
        cube = (x_index, y_index, z_index)

        # If the cube is not already in the dictionary, create a new group for it
        if cube not in groups:
            groups[cube] = []

        # append the atom to its respective cube
        row_dict = row.to_dict()
        groups[cube].append(row_dict)

    return groups, x_max, y_max, z_max
   # Print the groups of points
    # for i, (cube, group) in enumerate(groups.items()):
    #     print(f"Group {i+1} (Cube: {cube}):")
    #     for point in group:
    #         print(f"    {point}")


def adjacent_cubes(target, x_max, y_max, z_max):
    # Find the neighbors of a given cube
    neighbors = []

    x, y, z = target
    #print(f'Target={target}, x ={x}, y={y}, z={z}')
    for i in range(x-1, x+2):
        if i < 0: i = x_max
        if i > x_max: i = 0
        for j in range(y-1, y+2):
            if j < 0: j = y_max
            if j > y_max: j = 0
            for k in range(z-1, z+2):
                if k < 0: k = z_max
                if k > z_max: k = 0
                if (i, j, k) != target:
                    neighbors.append((i, j, k))
    #print(f'Neighbors: {neighbors}')
    return neighbors

# def get_neighbors(data):
#     attributes = ['ID', 'Type', 'X', 'Y', 'Z', 'Q', 'Neighbors', 'Distances']
#     neighbor_type = [1]
#     cutoff = 2.75
#     o_atoms = data.index[data['Type'].isin(neighbor_type)].tolist()

#     for index in o_atoms:
#         neighbors =[]
#         distances = []
#         target_atom = dict.fromkeys(attributes)

#         for col in attributes[:6]:
#              target_atom[f'{col}'] = data.loc[index][f'{col}']

#     return True

def build_neighbor_list(path, job, step, cutoff=2.5):
    # Splits the simulation into cubes and then builds a neighbor list by only searching current and adjacent cubes, returns a dictionary with additional 'Neighbors' and 'Distances' entries
    time, data, cell = ts.read_TS(path.format(job, step), return_cell=True)
            
    simulation, xm, ym, zm = group_atoms(data, cell)

    for cube in simulation:
        search_area = adjacent_cubes(cube, xm, ym, zm)
        search_area.append(cube)
        #print(f'For cube {cube}, the search area is {search_area}')
        for atom in simulation[cube]:
            atom['Neighbors'] = []
            atom['Distances'] = []
            target_xyz = [atom['X'], atom['Y'], atom['Z']]
        
            for area in search_area:
                for test_atom in simulation[area]:
                    test_xyz = [test_atom['X'], test_atom['Y'], test_atom['Z']]
                    bond_length = pbc_distance(target_xyz, test_xyz, cell)
                    if 0.005 < bond_length <= cutoff:
                        atom['Neighbors'].append(test_atom['Type'])
                        atom['Distances'].append(bond_length)
    return simulation


def pickle_simulation(data, pickle):
    # Takes the dictionary of cubes built from the simulation and generates a dataframe of just the atomic information and saves it as a pickle
    atom_list = []
    for cube in data:
        for atom in data[cube]:
            atom_list.append(atom)
   
    out_data = pd.DataFrame.from_records(atom_list)
    out_data.to_pickle(pickle)


#####################
def main():
    identifier = 'ni'
    rootdir = '/Users/steve/doe-project/recomb/ni_production'
    jobs = []

    # search rootdir for subdirectories and save them to a list
    for job in os.scandir(rootdir):
        if job.is_dir():
            name = str(job.path).split('/')[-1]
            if identifier in name:
                jobs.append(name)

    
    steps = range(0, 400001, 50000)
    cut = 3
    for job in tqdm(jobs):
        for step in steps:
            path = f'{rootdir}/{job}/chg.{step}.lmp'
            simulation = build_neighbor_list(path, job, step, cutoff=cut)
            print(f'For simulation {job} at time {step}')
            out_dir = '/Users/steve/doe-project/recomb/ni_production/pickles'
            pickle_simulation(simulation, f'{out_dir}/{job}/{job}_{step}.df')
            
if __name__ == "__main__":
    main()
