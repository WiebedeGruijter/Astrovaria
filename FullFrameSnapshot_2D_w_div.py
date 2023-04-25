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

def plot_2D_frame(dir, plot_dir, start_nr, end_nr, xy=False, r_p=False, zoom_level=None, vmin=None, vmax=None):

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

    list = list[start_nr:end_nr]
    print(list)

    n = 0 + start_nr
    for i in range(len(list)):

        fig = plt.figure(figsize=(5, 5))
        ax = plt.subplot(111)
        lim = 2e12
        
        snapshot = list[i]

        # read snapshot
        dblank = pw.ar.athdf(snapshot, quantities=[], level=mylevel, subsample=True)

        if xy == True:
            #x1_max = r[cm], x2_max = theta[rad], x3_max = phi[rad]
            x2sliceval = dblank['x2v'][np.argmin(np.abs(dblank['x2v'] - np.pi / 2))]
            d = pw.read_data(snapshot, orb, level=mylevel, x2_max=x2sliceval, x2_min=x2sliceval)
            r = pw.get_plot_array_midplane(d['gx1v'][:, 0, :])
            th = pw.get_plot_array_midplane(d['gx2v'][:, 0, :])
            phi = pw.get_plot_array_midplane(d['gx3v'][:, 0, :])
            v_r = pw.get_plot_array_midplane(d['vel1'][:, 0, :]) #Radial and azimuthal coordinates defined on a spherical grid.
            v_phi = pw.get_plot_array_midplane(d['vel3'][:, 0, :])

            X = r*np.sin(th)*np.cos(phi) / lenscale
            Y = r*np.sin(th)*np.sin(phi) /lenscale
            Z = r**-2 * np.gradient(r**2 * v_r, axis=1)/np.gradient(r, axis=1) + 1/(r*np.sin(th)) * np.gradient(v_phi, axis=0)/np.gradient(phi, axis=0)

            filename = 'div'+'_xy'
        else:
            # gives only left side of the full frame in xz plane
            x3sliceval = dblank['x3v'][np.argmin(np.abs(dblank['x3v']- np.pi))]
            d = pw.read_data(snapshot, orb, level=mylevel, x3_min=x3sliceval, x3_max=x3sliceval)
            X = pw.get_plot_array_midplane(d['x'][0, :, :]) / lenscale
            Y = pw.get_plot_array_midplane(d['z'][0, :, :]) / lenscale
            Z = pw.get_plot_array_midplane(np.log10(d[dvar][0, :, :] * SCALE))
            filename = 'div'+'_xz'

        if zoom_level!=None:
            im = ax.pcolormesh(X+a/lenscale, Y, Z, cmap='magma', vmin=vmin, vmax=vmax, shading='nearest', rasterized=True)
            lim = (zoom_level*lenscale)
            filename = filename+'_zoomed'
        else:
            im = ax.pcolormesh(X, Y, Z, cmap='magma', vmin=vmin, vmax=vmax, shading='nearest', rasterized=True)

        ax.set_aspect('equal')
        ax.set_xlim(-lim / lenscale, lim / lenscale)
        ax.set_ylim(-lim / lenscale, lim / lenscale)

        # sub region of the original image
        if zoom_level==None:
            axins = ax.inset_axes([0.57, 0.02, 0.43, 0.4])
            axins.pcolormesh(X,Y,Z, cmap='magma', vmin=vmin, vmax=vmax, shading='nearest', rasterized=True)

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

        cax = fig.add_axes([0.96, 0.21, 0.01, 0.6])
        
        cb = fig.colorbar(im, cax=cax, label=r'$\log_{10}\left(div \right) \ \ [{\rm g \ cm}^{-3}]$', extend='both')
            
        if n<10:
            nr_zeros = '000'
        elif (n>9) and (n<100):
            nr_zeros='00'
        elif (n>99) and (n<1000):
            nr_zeros='0'

        subfolder = ''
        if zoom_level!=None:
            subfolder='zoomed/'

        plt.savefig(plot_dir + subfolder + nr_zeros + str(n) + '_FullFrameSnapshot_'+filename+'.png', bbox_inches='tight', dpi=300)
        #plt.show()
        plt.close()

        print('plot ' + str(n) + ' done')
        n += 1

    print('All output files from ' + str(start_nr) + ' to ' + str(n - 1) + ' have been plotted.' )

#######################################################################################################################

# define input and parameters
f = '0.05_DN_2D_V8_l5_100_e0.20'                                   # name of folder
    
dir = '/Users/wiebe/Documents/Running-Athena/' + f + '/'       # path to folder
plot_dir = dir + 'plots/div/'                                   # save plots here

a = 7.47e11      # orbital distance in cm -- a = 0.05 au
rp = 6.71e9        # planetary radius in cm -- WASP-107b

SCALE = 1          # global scaling parameter for pressure and radius

appendix = 'pws_3D'
orb = pw.read_trackfile(dir + 'pm_trackfile.dat')

start_nr = 255        # selection of snapshot files in the directory, if only one snapshot: 0, 1
end_nr = 256

mylevel = 3         # resolution level, choose between 1 and 5


plot_2D_frame(dir, plot_dir, start_nr, end_nr, xy=True, r_p=True, zoom_level=None, vmin=None, vmax=None)

