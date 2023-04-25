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

def plot_2D_frame(dir, plot_dir, start_nr, xy=False, r_p=False, dvar='rho', zoom_level=None, vmin=None, vmax=None):

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

    list = list[start_nr:]
    print(list)

    snapshot = list[0]
    # read snapshot
    dblank = pw.ar.athdf(snapshot, quantities=[], level=mylevel, subsample=True)

    x2sliceval = dblank['x2v'][np.argmin(np.abs(dblank['x2v'] - np.pi / 2))]
    d = pw.read_data(snapshot, orb, level=mylevel, x2_max=x2sliceval, x2_min=x2sliceval)

    xp, yp, zp = pw.pos_secondary(orb, d['Time'])

    dp = np.sqrt(xp**2+yp**2)

    d2 = -1
    i=1

    while d2<dp: #Find the closest point of planet to star
        
        if d2>0:
            dp = d2

        snapshot = list[i]

        # read snapshot
        dblank = pw.ar.athdf(snapshot, quantities=[], level=mylevel, subsample=True)

        x2sliceval = dblank['x2v'][np.argmin(np.abs(dblank['x2v'] - np.pi / 2))]
        d = pw.read_data(snapshot, orb, level=mylevel, x2_max=x2sliceval, x2_min=x2sliceval)

        x2, y2, z2 = pw.pos_secondary(orb, d['Time'])
        d2 = np.sqrt(x2**2+y2**2)
        
        print('Distance planet-star:', dp)
        i+=1


    X = pw.get_plot_array_midplane(d['x'][:, 0, :]) 
    Y = pw.get_plot_array_midplane(d['y'][:, 0, :]) 
    Z = pw.get_plot_array_midplane((d[dvar][:, 0, :] * SCALE))

    # while d2>dp: #Find the farthest point of planet to star

    #     if d2>0:
    #         dp = d2

    #     snapshot = list[i]

    #     # read snapshot
    #     dblank = pw.ar.athdf(snapshot, quantities=[], level=mylevel, subsample=True)

    #     x2sliceval = dblank['x2v'][np.argmin(np.abs(dblank['x2v'] - np.pi / 2))]
    #     d = pw.read_data(snapshot, orb, level=mylevel, x2_max=x2sliceval, x2_min=x2sliceval)

    #     x2, y2, z2 = pw.pos_secondary(orb, d['Time'])
    #     d2 = np.sqrt(x2**2+y2**2)
        
    #     print('here:', dp)
    #     i+=1

    
    # read snapshot
    dblank = pw.ar.athdf(snapshot, quantities=[], level=mylevel, subsample=True)

    x2sliceval = dblank['x2v'][np.argmin(np.abs(dblank['x2v'] - np.pi / 2))]
    d = pw.read_data(snapshot, orb, level=mylevel, x2_max=x2sliceval, x2_min=x2sliceval)
    X = pw.get_plot_array_midplane(d['x'][:, 0, :]) / lenscale
    Y = pw.get_plot_array_midplane(d['y'][:, 0, :]) / lenscale
    Z = pw.get_plot_array_midplane((d[dvar][:, 0, :] * SCALE))

    
    fig = plt.figure(figsize=(5, 5))
    ax = plt.subplot(111)
    lim = 2e12
    filename = dvar

        
    if zoom_level!=None:
        im = ax.pcolormesh(X+a/lenscale, Y, Z, cmap='bwr', vmin=vmin, vmax=vmax, shading='nearest', rasterized=True)
        lim = (zoom_level*lenscale)
        filename = filename+'_zoomed'
    else:
        im = ax.pcolormesh(X, Y, Z, cmap='bwr', vmin=vmin, vmax=vmax, rasterized=True)

    ax.set_aspect('equal')
    ax.set_xlim(-lim / lenscale, lim / lenscale)
    ax.set_ylim(-lim / lenscale, lim / lenscale)

    # sub region of the original image
    if zoom_level==None:
        axins = ax.inset_axes([0.57, 0.02, 0.43, 0.4])
        axins.pcolormesh(X,Y,Z, cmap='magma', vmin=vmin, vmax=vmax, rasterized=True)

        ilim = 3e11               # size of subregion
        axins.set_aspect('equal')
        axins.set_xlim((-ilim - a) / lenscale, (ilim - a) / lenscale)
        axins.set_ylim(-ilim / lenscale, ilim / lenscale)
        axins.set_xticklabels('')
        axins.set_yticklabels('')
        ax.indicate_inset_zoom(axins)

    # labels and colorbar
    ax.set_xlabel(r'$x$ '+unit_name)
    if xy == True:
        ax.set_ylabel(r'$y$ '+unit_name)
    else: ax.set_ylabel(r'$z$ '+unit_name)

    # place a text box in upper left in axes coords
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.05, 0.95, 't = '+str(np.round(d['Time']/24/3600, 2))+' days', transform=ax.transAxes, fontsize=14,
    verticalalignment='top', bbox=props)

    #Place the planet in the frame
    # x2, y2, z2 = pw.pos_secondary(orb, d['Time'])
    ax.scatter(x2/lenscale, y2/lenscale, marker='o', s=10, color='green')

    cax = fig.add_axes([0.96, 0.21, 0.01, 0.6])
    if dvar == 'rho':
        cb = fig.colorbar(im, cax=cax, label=r'$\log_{10}\left(\rho \right) \ \ [{\rm g \ cm}^{-3}]$', extend='both')
        
    else:
        cb = fig.colorbar(im, cax=cax, label=r'$\log_{10}\left({\rm P} \right) \ \ [{\rm barye}]$', extend='both')

    subfolder = ''
    if zoom_level!=None:
        subfolder='zoomed/'


    plt.savefig(plot_dir + dvar +'/' + subfolder + '_FullFrameSnapshot_'+filename+'.png', bbox_inches='tight', dpi=300)
    #plt.show()
    plt.close()

    print('plot done:', plot_dir + dvar +'/' + subfolder + '_FullFrameSnapshot_'+filename)


#######################################################################################################################

# define input and parameters
f = '0.05_DN_2D_V8_l5_100_e0.20'                                    # name of folder
dir = '/Users/wiebe/Documents/Running-Athena/' + f + '/'       # path to folder
plot_dir = dir + 'plots/report/'                                   # save plots here

a = 7.47e11      # orbital distance in cm -- a = 0.05 au
rp = 6.71e9        # planetary radius in cm -- WASP-107b

SCALE = 1          # global scaling parameter for pressure and radius

appendix = 'pws_3D'
orb = pw.read_trackfile(dir + 'pm_trackfile.dat')




start_nr = 63        # selection of snapshot files in the directory, if only one snapshot: 0, 1


mylevel = 3         # resolution level, choose between 1 and 5

# dict_keys(['Coordinates', 'DatasetNames', 'MaxLevel', 'MeshBlockSize', 'NumCycles', 'NumMeshBlocks', 'NumVariables', 
# 'RootGridSize', 'RootGridX1', 'RootGridX2', 'RootGridX3', 'Time', 'VariableNames', 
# 'x1f', 'x1v', 'x2f', 'x2v', 'x3f', 'x3v', 'rho', 'press', 'vel1', 'vel2', 'vel3', 
# 'r0', 'r1', 'gx1v', 'gx2v', 'gx3v', 'dA', 'dvol', 'x', 'y', 'z', 'vx', 'vy', 'vz'])


plot_2D_frame(dir, plot_dir, start_nr, xy=True, r_p=True, dvar='vx', zoom_level=None, vmin=-2.5E7, vmax=2.5E7)

