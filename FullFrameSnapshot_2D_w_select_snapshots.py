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
plt.rcParams['font.size'] = 24

#######################################################################################################################

def plot_2D_frame(dir, plot_dir, xy=False, r_p=False, dvar='rho', zoom_level=None, vmin=None, vmax=None):

    if r_p is False:
        lenscale = c.au
        unit_name = '[au]'
    else:
        lenscale = rp
        unit_name = r'$[r_p]$'

    # list of all '.athdf' files in directory
    list = []
    for file_dir in glob2.glob(dir + '*.athdf'):
        list.append(file_dir)
    list.sort()

    list = [list[i] for i in numbers]
    print(list)

    fig, ax = plt.subplots(1, len(list), figsize=((7*len(list), 7)), sharey=True)
    for i in range(len(list)):

        lim = 2e12
        filename = dvar
        
        snapshot = list[i]

        # read snapshot
        dblank = pw.ar.athdf(snapshot, quantities=[], level=mylevel, subsample=True)

        if xy == True:
            #x1_max = r[cm], x2_max = theta[rad], x3_max = phi[rad]
            x2sliceval = dblank['x2v'][np.argmin(np.abs(dblank['x2v'] - np.pi / 2))]
            d = pw.read_data(snapshot, orb, level=mylevel, x2_max=x2sliceval, x2_min=x2sliceval)
            X = pw.get_plot_array_midplane(d['x'][:, 0, :]) / lenscale
            Y = pw.get_plot_array_midplane(d['y'][:, 0, :]) / lenscale
            if dvar=='rho':
                Z = pw.get_plot_array_midplane((np.log10(d[dvar][:, 0, :] * SCALE)))
            else:
                Z = pw.get_plot_array_midplane((d[dvar][:, 0, :] * SCALE))
            filename = dvar+'_xy'
        else:
            # gives only left side of the full frame in xz plane
            x3sliceval = dblank['x3v'][np.argmin(np.abs(dblank['x3v']- np.pi))]
            d = pw.read_data(snapshot, orb, level=mylevel, x3_min=x3sliceval, x3_max=x3sliceval)
            X = pw.get_plot_array_midplane(d['x'][0, :, :]) / lenscale
            Y = pw.get_plot_array_midplane(d['z'][0, :, :]) / lenscale
            Z = pw.get_plot_array_midplane((d[dvar][0, :, :] * SCALE))
            filename = dvar+'_xz'

        x2, y2, z2 = pw.pos_secondary(orb, d['Time']) 

        if zoom_level!=None:
            im = ax[i].pcolormesh(X-x2/lenscale, Y-y2/lenscale, Z, cmap='magma', vmin=vmin, vmax=vmax, shading='nearest', rasterized=True)
            lim = (zoom_level*lenscale)
            filename = filename+'_zoomed'
            ax[i].set_xticks([-zoom_level/2, 0, zoom_level/2])
            ax[i].set_yticks([-zoom_level/2, 0, zoom_level/2])
        else:
            im = ax[i].pcolormesh(X, Y, Z, cmap='magma', vmin=vmin, vmax=vmax, rasterized=True)

        ax[i].set_aspect('equal')
        ax[i].set_xlim(-lim / lenscale, lim / lenscale)
        ax[i].set_ylim(-lim / lenscale, lim / lenscale)

        # sub region of the original image
        if zoom_level==None:
            axins = ax[i].inset_axes([0.57, 0.02, 0.43, 0.4])
            axins.pcolormesh(X,Y,Z, cmap='bwr', vmin=vmin, vmax=vmax, rasterized=True)

            ilim = 3e11               # size of subregion
            axins.set_aspect('equal')
            axins.set_xlim((-ilim - a) / lenscale, (ilim - a) / lenscale)
            axins.set_ylim(-ilim / lenscale, ilim / lenscale)
            axins.set_xticklabels('')
            axins.set_yticklabels('')
            ax[i].indicate_inset_zoom(axins)

        # labels and colorbar
        ax[i].set_xlabel(r'$x$ '+unit_name)
        if xy == True:
            ax[0].set_ylabel(r'$y$ '+unit_name)
        else: ax[0].set_ylabel(r'$z$ '+unit_name)
        

        # place a text box in axes coords to display the phase of the orbit
        time = d['Time'] #All snapshots are taken at the same time, so we can simply take the time of the last one
        period = np.sqrt(4*np.pi**2*a**3/(6.67430e-8)/(Mp+Mstar))
        phase = time/period + 0.5 #We add 0.5 since we start our orbit at phi=0.5
        phase = np.round(phase, 2)
        while phase>=1:
            phase = phase-1
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax[i].text(0.05, 0.95, 'phase = '+str(np.round(phase, 2)), transform=ax[i].transAxes, fontsize=24,
        verticalalignment='top', bbox=props)


        #Plot a circle around the planet
        #x2, y2, z2 = pw.pos_secondary(orb, d['Time']) 
        if zoom_level==None:
            circle1 = plt.Circle((x2/lenscale, y2/lenscale), rp/lenscale, fill=False, color='green', linewidth=0.5)
            circle2 = plt.Circle((x2/lenscale, y2/lenscale), 10*rp/lenscale, fill=False, color='green', linewidth=0.5)
        else: 
            circle1 = plt.Circle((0,0), rp/lenscale, fill=False, color='green', linewidth=0.5)
            circle2 = plt.Circle((0,0), 10*rp/lenscale, fill=False, color='green', linewidth=1.0)
            circle3 = plt.Circle((0,0), 5*rp/lenscale, fill=False, color='green', linewidth=1.0)
        ax[i].add_patch(circle1)
        ax[i].add_patch(circle2)
        ax[i].add_patch(circle3)

        cax = fig.add_axes([0.92, 0.21, 0.01, 0.6])
        
        if dvar == 'rho':
            cb = fig.colorbar(im, cax=cax, label=r'$\log_{10}\left(\rho \right) \ \ [{\rm g \ cm}^{-3}]$', extend='both', ticks=[vmin, vmax, np.mean([vmin, vmax])])
            
        else:
            cb = fig.colorbar(im, cax=cax, label=r'$V_{los}$[cm/s]', extend='both')
            
        subfolder = ''
        if zoom_level!=None:
            subfolder='zoomed/'

    plt.subplots_adjust(wspace=0, hspace=0)
    plt.savefig(plot_dir + dvar +'/' + subfolder + '_FullFrameSnapshot_'+filename+f[-4:]+'.png', bbox_inches='tight', dpi=300)
    #plt.show()
    plt.close()

        
        

    print('All output files have been plotted.' )

#######################################################################################################################

# define input and parameters
f = '0.05_DN_2D_V8_l5_100_e0.20'                                    # name of folder
dir = '/Users/wiebe/Documents/Running-Athena/' + f + '/'       # path to folder
plot_dir = dir + 'plots/multi/'                                   # save plots here

a = 7.47e11      # orbital distance in cm -- a = 0.05 au
rp = 6.71e9        # planetary radius in cm -- WASP-107b
Mp = 1.82e29        # planetary mass in gram
Mstar = 1.36e33     # stellar mass in gram

SCALE = 1          # global scaling parameter for pressure and radius

appendix = 'pws_3D'
orb = pw.read_trackfile(dir + 'pm_trackfile.dat')




#numbers = [234, 245, 256, 266, 276]
numbers = [234, 245, 256, 266]

mylevel = 3         # resolution level, choose between 1 and 5


# dict_keys(['Coordinates', 'DatasetNames', 'MaxLevel', 'MeshBlockSize', 'NumCycles', 'NumMeshBlocks', 'NumVariables', 
# 'RootGridSize', 'RootGridX1', 'RootGridX2', 'RootGridX3', 'Time', 'VariableNames', 
# 'x1f', 'x1v', 'x2f', 'x2v', 'x3f', 'x3v', 'rho', 'press', 'vel1', 'vel2', 'vel3', 
# 'r0', 'r1', 'gx1v', 'gx2v', 'gx3v', 'dA', 'dvol', 'x', 'y', 'z', 'vx', 'vy', 'vz'])

        # if dvar=='vel1':
        #     colormap = 'bwr_r'

plot_2D_frame(dir, plot_dir, xy=True, r_p=True, dvar='rho', zoom_level=40, vmin=-22, vmax=-16)


