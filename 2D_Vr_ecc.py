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

def plot_2D_frame(dir, plot_dir, start_nr, end_nr, xy=False, r_p=False, dvar='rho', zoom_level=None, vmin=None, vmax=None):

    if r_p is False:
        lenscale = c.au
        unit_name = '[au]'
    else:
        lenscale = rp
        unit_name = r'$[r_p]$'

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

    n = 0 + start_nr
    for i in range(len(list)):

        fig, ax = plt.subplots(1, 1, figsize=((15,15)), sharey=True)
        lim = 2e12
        filename = dvar

        vel1_list = []
        
        for j in range(len(tot_list)):
            snapshot = tot_list[j][i]

            orb = pw.read_trackfile(dir+f[j]+'/' + 'pm_trackfile.dat')

            # read snapshot
            dblank = pw.ar.athdf(snapshot, quantities=[], level=mylevel, subsample=True)

            
            #x1_max = r[cm], x2_max = theta[rad], x3_max = phi[rad]
            x2sliceval = dblank['x2v'][np.argmin(np.abs(dblank['x2v'] - np.pi / 2))]
            d = pw.read_data(snapshot, orb, level=mylevel, x2_max=x2sliceval, x2_min=x2sliceval)
            X = pw.get_plot_array_midplane(d['x'][:, 0, :]) 
            Y = pw.get_plot_array_midplane(d['y'][:, 0, :]) 
            Z = pw.get_plot_array_midplane((d[dvar][:, 0, :] * SCALE))
            
            #Position of the planet
            x2, y2, z2 = pw.pos_secondary(orb, d['Time']) 

            mean_vr = np.mean(Z[(((X-x2)**2+(Y-y2)**2)<((10*lenscale)**2))&(((X-x2)**2+(Y-y2)**2)>((1*lenscale)**2))])
            # mean_vr = np.mean(Z[((X-x2)**2+(Y-y2)**2)<((10*lenscale)**2)])
            
            vel1_list.append(mean_vr)
        
        
            
        # place a text box in upper left in axes coords to display the phase of the orbit
        time = d['Time'] #All snapshots are taken at the same time, so we can simply take the time of the last one
        period = np.sqrt(4*np.pi**2*a**3/(6.67430e-8)/(Mp+Mstar))
        phase = time/period + 0.5 #We add 0.5 since we start our orbit at phi=0.5
        while phase>1:
            phase = phase-1
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(0.5, 0.95, 'phase = '+str(np.round(phase, 2)), transform=ax.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)

        plt.plot(eccentricities, vel1_list)
        plt.scatter(eccentricities, vel1_list, color='orange')
        plt.xticks(eccentricities)
        #plt.ylim(-3e6, 3e6)
        plt.ylim(0, 4e-18)
        plt.title('Radial velocity averaged 10 $R_p$ around planet')
        plt.ylabel('$V_r$ [cm/s]')
        plt.xlabel('Eccentricity')

        if n<10:
            nr_zeros = '000'
        elif (n>9) and (n<100):
            nr_zeros='00'
        elif (n>99) and (n<1000):
            nr_zeros='0'
        plt.subplots_adjust(wspace=0)
        plt.savefig(plot_dir + dvar +'/'+ nr_zeros + str(n) + '_FullFrameSnapshot_'+filename+'.png', bbox_inches='tight', dpi=300)
        #plt.show()
        plt.close()

        print('plot ' + str(n) + ' done')
        n += 1

    print('All output files from ' + str(start_nr) + ' to ' + str(n - 1) + ' have been plotted.' )

#######################################################################################################################

# define input and parameters
f = ['0.05_DN_2D_V8_I5_100_e0.00', '0.05_DN_2D_V8_I5_100_e0.05', '0.05_DN_2D_V8_I5_100_e0.10', '0.05_DN_2D_V8_I5_100_e0.15', '0.05_DN_2D_V8_l5_100_e0.20']                                    # name of folder
dir = '/Users/wiebe/Documents/Running-Athena/'    # path to folder

plot_dir = dir + 'Ecc_to_var/'                                   # save plots here

eccentricities = [0, 0.05, 0.10, 0.15, 0.20]

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


plot_2D_frame(dir, plot_dir, start_nr, end_nr, xy=True, r_p=True, dvar='rho', zoom_level=None, vmin=None, vmax=None)

