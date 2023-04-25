from planet_wind_constants import c
import planet_wind_utils_v6 as pw
import numpy as np
import matplotlib.pyplot as plt
import glob2

# set some global options for plots
plt.rcParams['figure.figsize'] = (6, 5)
plt.rcParams['legend.frameon'] = True
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['legend.borderpad'] = 0.2
plt.rcParams['legend.labelspacing'] = 0.2
plt.rcParams['legend.handletextpad'] = 0.2
plt.rcParams['font.family'] = 'stixgeneral'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.size'] = 16

#######################################################################################################################

def plot_2D_frame(dir, plot_dir, start_nr, end_nr, xy=False, r_p=False, dvar='rho', inner_circle=1, outer_circle=10):

    if r_p is False:
        lenscale = c.au
        unit_name = '[au]'
    else:
        lenscale = rp
        unit_name = r' $r_p$'

    #Here we make a list of lists that contains a list of the hdf files for each different folder
    tot_list = []
    # list of all '.athdf' files in directory
    for filename in f:
   
        list = []
        for file_dir in glob2.glob(dir+filename+'/' + '*.athdf'):
            list.append(file_dir)
        list.sort()

        list = list[start_nr:end_nr]
        tot_list.append(list)

    print(tot_list)

    fig, ax = plt.subplots(1, 1, figsize=((15,15)), sharey=True)
    lim = 2e12
    filename = dvar

    
    phase_axis = np.zeros(len(list))
    vel1_axis = np.zeros((len(tot_list), len(list)))
    r_axis = np.zeros((len(tot_list), len(list)))
    time_axis = np.zeros((len(tot_list), len(list)))
    for i in range(len(list)):
        for j in range(len(tot_list)):
            print(j, i)
            snapshot = tot_list[j][i]

            orb = pw.read_trackfile(dir+f[j]+'/' + 'pm_trackfile.dat')

            # read snapshot
            dblank = pw.ar.athdf(snapshot, quantities=[], level=mylevel, subsample=True)

            
            #x1_max = r[cm], x2_max = theta[rad], x3_max = phi[rad]
            x2sliceval = dblank['x2v'][np.argmin(np.abs(dblank['x2v'] - np.pi / 2))]
            d = pw.read_data(snapshot, orb, level=mylevel, x2_max=x2sliceval, x2_min=x2sliceval)
            X = pw.get_plot_array_midplane(d['x'][:, 0, :]) 
            Y = pw.get_plot_array_midplane(d['y'][:, 0, :]) 
            rho = pw.get_plot_array_midplane((d['rho'][:, 0, :] * SCALE))
            press = pw.get_plot_array_midplane((d['press'][:, 0, :] * SCALE))
            
            gamma = 1.0001
            mu = 0.6
            Z = gamma * mu * c.mp * press / (rho * c.kB)
            
            #Position of the planet
            x2, y2, z2 = pw.pos_secondary(orb, d['Time']) 

            mean_vr = np.mean(Z[(((X-x2)**2+(Y-y2)**2)<((outer_circle*lenscale)**2))&(((X-x2)**2+(Y-y2)**2)>((inner_circle*lenscale)**2))])
            
            vel1_axis[j, i] = mean_vr

         
        # calculate phase
        time = d['Time'] # All snapshots are taken at the same time, so we can simply take the time of the last one
        period = np.sqrt(4*np.pi**2*a**3/(6.67430e-8)/(Mp+Mstar)) # Period from Kepler's third law
        phase = time/period + 0.5 # We add 0.5 since we start our orbit at phi=0.5
        phase_axis[i] = phase

        print('File '+str(i)+' read')
    
    while phase_axis[-1]>1: #Shift the phase so that the final plotted orbit ranges from phase=0 to phase=1
        phase_axis+=-1
    if phase_axis[0]<0:
        phase_axis+=1

    for j in range(len(tot_list)):
        
        plt.plot(phase_axis, vel1_axis[j], label='T for e = '+str(eccentricities[j]), color=colors[j])

        #Plot the radial velocity of the planet
        # vr_planet = np.gradient(r_axis[j])/np.gradient(time_axis[j])
        # plt.plot(phase_axis[1:-1], vr_planet[1:-1], label='Planet velocity for e = '+str(eccentricities[j]), linestyle='dotted', color=colors[j])

    
    #plt.ylim(-3e6, 7e6)
    plt.title('Radial velocity averaged from '+str(inner_circle)+unit_name+' to '+str(outer_circle)+unit_name+ ' around centre of planet')
    plt.ylabel('$V_r$ [cm/s]')
    plt.xlabel('Phase')
    plt.legend()
    plt.subplots_adjust(wspace=0)
    plt.savefig(plot_dir + dvar +'/'+ '_FullFrameSnapshot_'+filename+'_'+str(inner_circle)+'to'+str(outer_circle)+'.png', bbox_inches='tight', dpi=300)
    #plt.show()
    plt.close()

    print('All output files have been plotted.' )

#######################################################################################################################

# define input and parameters
f = ['0.05_DN_2D_V8_I5_100_e0.00', '0.05_DN_2D_V8_I5_100_e0.05', '0.05_DN_2D_V8_I5_100_e0.10', '0.05_DN_2D_V8_I5_100_e0.15', '0.05_DN_2D_V8_l5_100_e0.20']                                    # name of folder
dir = '/Users/wiebe/Documents/Running-Athena/'    # path to folder

plot_dir = dir + 'Phi_to_var/'                                   # save plots here

colors = ['b', 'g', 'r', 'c', 'm']

eccentricities = [0, 0.05, 0.10, 0.15, 0.20]    # list of eccentricities corresponding to the directories in f


a = 7.47e11      # orbital distance in cm -- a = 0.05 au
rp = 6.71e9        # planetary radius in cm -- WASP-107b
Mp = 1.82e29        # planetary mass in gram
Mstar = 1.36e33     # stellar mass in gram

SCALE = 1          # global scaling parameter for pressure and radius

appendix = 'pws_3D'

start_nr = 235        # selection of snapshot files in the directory, if only one snapshot: 0, 1
end_nr = 276

mylevel = 3         # resolution level, choose between 1 and 5


# dict_keys(['Coordinates', 'DatasetNames', 'MaxLevel', 'MeshBlockSize', 'NumCycles', 'NumMeshBlocks', 'NumVariables', 
# 'RootGridSize', 'RootGridX1', 'RootGridX2', 'RootGridX3', 'Time', 'VariableNames', 
# 'x1f', 'x1v', 'x2f', 'x2v', 'x3f', 'x3v', 'rho', 'press', 'vel1', 'vel2', 'vel3', 
# 'r0', 'r1', 'gx1v', 'gx2v', 'gx3v', 'dA', 'dvol', 'x', 'y', 'z', 'vx', 'vy', 'vz'])


plot_2D_frame(dir, plot_dir, start_nr, end_nr, xy=True, r_p=True, dvar='T', inner_circle=1, outer_circle=10)

