# Plot the displacement of the cathode and anode working atoms:
import math, os
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np

def read_TS(file): # read in the data file for a given timestep and store it as a list of dictionaries, each dict represents an atom
    input = open(file, 'r')
    raw_data = input.readlines()
    input.close()

    timestep = int(raw_data[1])
   
    run_data = []
    for line in raw_data[9:]:
                ld = line.split()
                ts_data = dict.fromkeys(['ID', 'Type', 'X', 'Y', 'Z', 'Q'])
                #print(ld)
                ts_data['ID'] = int(ld[0])
                ts_data['Type'] = int(ld[1])
                ts_data['X'] = float(ld[2])
                ts_data['Y'] = float(ld[3])
                ts_data['Z'] = float(ld[4])
                ts_data['Q'] = float(ld[5])
                run_data.append(ts_data)
    run_data = sorted(run_data, key=lambda x: x['ID']) # sort by atom ID tag

    cell_x = float(raw_data[5].split()[0]) + float(raw_data[5].split()[1]) # casts from string in scientific notation to float and sums lammps xmin & xmax
    cell_y = float(raw_data[6].split()[0]) + float(raw_data[6].split()[1])
    cell_z = float(raw_data[7].split()[0]) + float(raw_data[7].split()[1])
    cell = [cell_x, cell_y, cell_z]
    return timestep, run_data, cell

def pbc_distance(a, b, cell):# Calculate the distance between atom a(list with [x,y,z]) and atom b in a periodic cell with dimensions cell=[x,y,z]
    A = cell[0]
    B = cell[1]
    C = cell[2]

    dx = abs(a[0] - b[0])
    x = min(dx, abs(A - dx))
     
    dy = abs(a[1] - b[1])
    y = min(dy, abs(B - dy))
     
    dz = abs(a[2] - b[2])
    z = min(dz, abs(C - dz))
 
    return math.sqrt(x**2 + y**2 + z**2)

def msd(type, path, t1, t0): # return the root mean squared deviation of atoms of type given two timesteps
    time1, data1, cell1 = read_TS(f'{path}/chg.{t1}.lmp')
    data1 = [atom for atom in data1 if atom['Type'] == type] # only save atoms of the desired type for t1

    time0, data0, cell0 = read_TS(f'{path}/chg.{t0}.lmp') 
    data0 = [atom for atom in data0 if atom['Type'] == type]
    
    # iterate over every atom in data0 by index so we can access the same atom in data1
    msd = 0
    for i in range(0, len(data0)):
        #print(f"{data0[i]}, {data1[i]}")
        r1 = [data1[i]['X'], data1[i]['Y'], data1[i]['Z']]
        r0 = [data0[i]['X'], data0[i]['Y'], data0[i]['Z']]
        dist = pbc_distance(r1, r0, cell1)
        msd += dist**2
    msd = math.sqrt(msd / len(data0))
    
    return msd

def get_data(root_dir, job, cathode_id, anode_id):
    timesteps = []
    cat_disp = []
    an_disp = []
    for file in os.listdir(f'{root_dir}/{job}'):
        if file.startswith('chg'):
            ts = int(file.split('.')[1])
            timesteps.append(ts)

    timesteps.sort()
    print(f'Plotting {job}\n')
    for t in tqdm(range(1, len(timesteps))):
        m1 = msd(cathode_id, f'{root_dir}/{job}', timesteps[t], timesteps[0])
        cat_disp.append(m1)
        m2 = msd(anode_id, f'{root_dir}/{job}', timesteps[t], timesteps[0])
        an_disp.append(m2)
    return cat_disp, an_disp, 

def plot_slab(cathode, anode):
    timesteps = [x*2.5 for x in range(0,201,1) ] #output data every 10,000 timesteps = every 2.5ps for 500ps
    fig, axs = plt.subplots(3,2)
    
    # Plot Pt RMSD Plots
        # Plot Pt Cathodes
    axs[0,0].plot(timesteps[1::25],cathode['pt']['0v'][1::25], marker='s', linestyle='dotted', color='black', label='0V')
    axs[0,0].plot(timesteps[5::25],cathode['pt']['1v'][5::25], marker='*', linestyle='dashed', color='black', label='1V') 
    axs[0,0].plot(timesteps[10::25],cathode['pt']['5v'][10::25], marker='+', linestyle='dashdot', color='black', label='5V')  
        # Plot Pt Anode
    axs[0,1].plot(timesteps[1::25],anode['pt']['0v'][1::25], marker='s', linestyle='dotted', color='black')
    axs[0,1].plot(timesteps[5::25],anode['pt']['1v'][5::25], marker='*', linestyle='dashed', color='black') 
    axs[0,1].plot(timesteps[10::25],anode['pt']['5v'][10::25], marker='+', linestyle='dashdot', color='black')
        # Plot Ni Cathodes
    trim_0v = timesteps[:len(cathode['ni']['0v'])]
    axs[1,0].plot(trim_0v[1::25],cathode['ni']['0v'][1::25], marker='s', linestyle='dotted', color='black')
    trim_1v = timesteps[:len(cathode['ni']['1v'])]
    axs[1,0].plot(trim_1v[5::25],cathode['ni']['1v'][5::25], marker='*', linestyle='dashed', color='black')
    trim_5v = timesteps[:len(cathode['ni']['5v'])] 
    axs[1,0].plot(trim_5v[10::25],cathode['ni']['5v'][10::25], marker='+', linestyle='dashdot', color='black')   
        # Plot Ni Cathodes
    trim_0v = timesteps[:len(anode['ni']['0v'])]
    axs[1,1].plot(trim_0v[1::25],anode['ni']['0v'][1::25], marker='s', linestyle='dotted', color='black')
    trim_1v = timesteps[:len(anode['ni']['1v'])]
    axs[1,1].plot(trim_1v[5::25],anode['ni']['1v'][5::25], marker='*', linestyle='dashed', color='black')
    trim_5v = timesteps[:len(anode['ni']['5v'])] 
    axs[1,1].plot(trim_5v[10::25],anode['ni']['5v'][10::25], marker='+', linestyle='dashdot', color='black')   
        # Plot Cu Cathodes
    trim_0v = timesteps[:len(cathode['cu']['0v'])]
    axs[2,0].plot(trim_0v[1::25],cathode['cu']['0v'][1::25], marker='s', linestyle='dotted', color='black')
    trim_1v = timesteps[:len(cathode['cu']['1v'])]
    axs[2,0].plot(trim_1v[5::25],cathode['cu']['1v'][5::25], marker='*', linestyle='dashed', color='black')
    trim_5v = timesteps[:len(cathode['cu']['5v'])] 
    axs[2,0].plot(trim_5v[10::25],cathode['cu']['5v'][10::25], marker='+', linestyle='dashdot', color='black')   
        # Plot Cu Cathodes
    trim_0v = timesteps[:len(anode['cu']['0v'])]
    axs[2,1].plot(trim_0v[1::25],anode['cu']['0v'][1::25], marker='s', linestyle='dotted', color='black')
    trim_1v = timesteps[:len(anode['cu']['1v'])]
    axs[2,1].plot(trim_1v[5::25],anode['cu']['1v'][5::25], marker='*', linestyle='dashed', color='black')
    trim_5v = timesteps[:len(anode['cu']['5v'])] 
    axs[2,1].plot(trim_5v[10::25],anode['cu']['5v'][10::25], marker='+', linestyle='dashdot', color='black')
    
    y_bound = (0.15, 0.65)
    plt.setp(axs, ylim=y_bound)
    axs[0,0].legend()

    axs[0,0].set_title('Cathodes')
    axs[0,1].set_title('Anodes')
    axs[2,0].set_xlabel('Simulation Time [ps]')
    axs[2,1].set_xlabel('Simulation Time [ps]')
    axs[0,0].set_ylabel('Root Mean Square Deviation')
    axs[1,0].set_ylabel('Root Mean Square Deviation')
    axs[2,0].set_ylabel('Root Mean Square Deviation')

    axs[0,1].text(500,0.4, 'Pt (111)')
    axs[1,1].text(500,0.4, 'Ni (111)')
    axs[2,1].text(500,0.4, 'Cu (111)')

    plt.show()
    # plt.savefig('pt_RMSD.png')
    plt.clf()


def plot_defect(cathode, anode):
    timesteps = [x*2.5 for x in range(0,203,1) ] #output data every 10,000 timesteps = every 2.5ps for 500ps
    fig, axs = plt.subplots(3,2)
    
    # Plot Pt RMSD Plots
        # Plot Pt Cathodes
    trim_0v = timesteps[:len(cathode['pt']['step'])]
    axs[0,0].plot(trim_0v[1::25],cathode['pt']['step'][::25], marker='s', linestyle='dotted', color='black', label='step-edge')
    trim_1v = timesteps[:len(cathode['pt']['ad'])]
    axs[0,0].plot(trim_1v[5::25],cathode['pt']['ad'][5::25], marker='*', linestyle='dashed', color='black', label='adatom') 
    trim_5v = timesteps[:len(cathode['pt']['hole'])]
    axs[0,0].plot(trim_5v[10::25],cathode['pt']['hole'][10::25], marker='+', linestyle='dashdot', color='black', label='hole')  
        # Plot Pt Anode
    trim_0v = timesteps[:len(anode['pt']['step'])]
    axs[0,1].plot(trim_0v[1::25],anode['pt']['step'][::25], marker='s', linestyle='dotted', color='black')
    trim_1v = timesteps[:len(anode['pt']['ad'])]
    axs[0,1].plot(trim_1v[5::25],anode['pt']['ad'][5::25], marker='*', linestyle='dashed', color='black') 
    trim_5v = timesteps[:len(anode['pt']['hole'])]
    axs[0,1].plot(trim_5v[10::25],anode['pt']['hole'][10::25], marker='+', linestyle='dashdot', color='black')
        # Plot Ni Cathodes
    trim_0v = timesteps[:len(cathode['ni']['step'])]
    axs[1,0].plot(trim_0v[1::25],cathode['ni']['step'][1::25], marker='s', linestyle='dotted', color='black')
    trim_1v = timesteps[:len(cathode['ni']['ad'])]
    axs[1,0].plot(trim_1v[5::25],cathode['ni']['ad'][5::25], marker='*', linestyle='dashed', color='black')
    trim_5v = timesteps[:len(cathode['ni']['hole'])] 
    axs[1,0].plot(trim_5v[10::25],cathode['ni']['hole'][10::25], marker='+', linestyle='dashdot', color='black')   
        # Plot Ni Anode
    trim_0v = timesteps[:len(anode['ni']['step'])]
    axs[1,1].plot(trim_0v[1::25],anode['ni']['step'][1::25], marker='s', linestyle='dotted', color='black')
    trim_1v = timesteps[:len(anode['ni']['ad'])]
    axs[1,1].plot(trim_1v[5::25],anode['ni']['ad'][5::25], marker='*', linestyle='dashed', color='black')
    trim_5v = timesteps[:len(anode['ni']['hole'])] 
    axs[1,1].plot(trim_5v[10::25],anode['ni']['hole'][10::25], marker='+', linestyle='dashdot', color='black')   
        # Plot Cu Cathodes
    trim_0v = timesteps[:len(cathode['cu']['step'])]
    axs[2,0].plot(trim_0v[1::25],cathode['cu']['step'][1::25], marker='s', linestyle='dotted', color='black')
    trim_1v = timesteps[:len(cathode['cu']['ad'])]
    axs[2,0].plot(trim_1v[5::25],cathode['cu']['ad'][5::25], marker='*', linestyle='dashed', color='black')
    trim_5v = timesteps[:len(cathode['cu']['hole'])] 
    axs[2,0].plot(trim_5v[10::25],cathode['cu']['hole'][10::25], marker='+', linestyle='dashdot', color='black')   
        # Plot Cu Anode
    trim_0v = timesteps[:len(anode['cu']['step'])]
    axs[2,1].plot(trim_0v[1::25],anode['cu']['step'][1::25], marker='s', linestyle='dotted', color='black')
    trim_1v = timesteps[:len(anode['cu']['ad'])]
    axs[2,1].plot(trim_1v[5::25],anode['cu']['ad'][5::25], marker='*', linestyle='dashed', color='black')
    trim_5v = timesteps[:len(anode['cu']['hole'])] 
    axs[2,1].plot(trim_5v[10::25],anode['cu']['hole'][10::25], marker='+', linestyle='dashdot', color='black')
    
    y_bound = (0.10, 1.25)
    plt.setp(axs, ylim=y_bound)
    axs[0,0].legend()

    axs[0,0].set_title('Cathodes')
    axs[0,1].set_title('Anodes')
    axs[2,0].set_xlabel('Simulation Time [ps]')
    axs[2,1].set_xlabel('Simulation Time [ps]')
    axs[0,0].set_ylabel('Root Mean Square Deviation')
    axs[1,0].set_ylabel('Root Mean Square Deviation')
    axs[2,0].set_ylabel('Root Mean Square Deviation')

    axs[0,1].text(575,0.4, 'Pt (111)')
    axs[1,1].text(575,0.4, 'Ni (111)')
    axs[2,1].text(575,0.4, 'Cu (111)')

    plt.show()
    # plt.savefig('pt_RMSD.png')
    plt.clf()
###################

def slab():
    plt_20 = False
    id = ['pt', 'ni', 'cu']
    types = ['separate']
    jobs = ['0v', '1v', '5v']
    cat_data = {key: None for key in id}
    anode_data = {key: None for key in id}
    for identifier in id:
        for jtype in types:
            cat_disp = {key: None for key in jobs}
            anode_disp = {key: None for key in jobs}
            for volt in jobs:
                root_dir = f'/Users/steve/doe-project/recomb/{identifier}_production'
                cathode_id = 4
                anode_id = 5
                work_dir = identifier+'_'+jtype+'_'+volt
                print(f'{root_dir}/{work_dir}')
                cat_disp[volt], anode_disp[volt] = get_data(root_dir, work_dir, cathode_id, anode_id)
        cat_data[identifier] = cat_disp
        anode_data[identifier] = anode_disp
    plot_slab(cat_data, anode_data)
            
def defect():
    plt_20 = False
    id = ['pt', 'ni', 'cu']
    jobs = ['step', 'ad', 'hole']
    cat_data = {key: None for key in id}
    anode_data = {key: None for key in id}
    for identifier in id:
            cat_disp = {key: None for key in jobs}
            anode_disp = {key: None for key in jobs}
            for volt in jobs:
                root_dir = f'/Users/steve/doe-project/recomb/defect_slabs/defect_fix'
                cathode_id = 4
                anode_id = 8
                work_dir = f'{identifier}-{volt}'
                print(f'{root_dir}/{work_dir}')
                cat_disp[volt], anode_disp[volt] = get_data(root_dir, work_dir, cathode_id, anode_id)

            cat_data[identifier] = cat_disp
            anode_data[identifier] = anode_disp
    plot_defect(cat_data, anode_data)





if __name__ == '__main__':
    slab()
    # defect()