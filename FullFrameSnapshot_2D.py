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

def plot_2D_frame(dir, plot_dir, start_nr, end_nr, xy=False):

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
            X = pw.get_plot_array_midplane(d['x'][:, 0, :]) / c.au
            Y = pw.get_plot_array_midplane(d['y'][:, 0, :]) / c.au
            RHO = pw.get_plot_array_midplane(np.log10(d['rho'][:, 0, :] * SCALE))
            #PRESS = pw.get_plot_array_midplane(np.log10(d['press'][:, 0, :] * SCALE))

        else:
            # gives only left side of the full frame in xz plane
            x3sliceval = dblank['x3v'][np.argmin(np.abs(dblank['x3v']- np.pi))]
            d = pw.read_data(snapshot, orb, level=mylevel, x3_min=x3sliceval, x3_max=x3sliceval)
            X = pw.get_plot_array_midplane(d['x'][0, :, :]) / c.au
            Y = pw.get_plot_array_midplane(d['z'][0, :, :]) / c.au
            RHO = pw.get_plot_array_midplane(np.log10(d['rho'][0, :, :] * SCALE))

        #im = ax.pcolormesh(X, Y, RHO, cmap='magma', vmin=-24, vmax=-16.5, shading='flat', rasterized=True)
        im = ax.pcolormesh(X, Y, RHO, cmap='magma', vmin=-24, vmax=-16.5, shading='nearest', rasterized=True)
        
        ax.set_aspect('equal')
        ax.set_xlim(-lim / c.au, lim / c.au)
        ax.set_ylim(-lim / c.au, lim / c.au)

        # sub region of the original image
        axins = ax.inset_axes([0.57, 0.02, 0.43, 0.4])
        #axins.pcolormesh(X,Y,RHO, cmap='magma', vmin=-24, vmax=-16.5, shading='flat', rasterized=True)
        axins.pcolormesh(X,Y,RHO, cmap='magma', vmin=-24, vmax=-16.5, shading='nearest', rasterized=True)

        ilim = 1.5e11               # size of subregion
        axins.set_aspect('equal')
        axins.set_xlim((-ilim - a) / c.au, (ilim - a) / c.au)
        axins.set_ylim(-ilim / c.au, ilim / c.au)
        axins.set_xticklabels('')
        axins.set_yticklabels('')
        ax.indicate_inset_zoom(axins)

        # labels and colorbar
        ax.set_xlabel(r'$x \ \ [\rm {au}]$')
        if xy == True:
            ax.set_ylabel(r'$y \ \ [\rm {au}]$')
        else: ax.set_ylabel(r'$z \ \ [\rm {au}]$')

        cax = fig.add_axes([0.96, 0.21, 0.01, 0.6])
        cb = fig.colorbar(im, cax=cax, label=r'$\log_{10}\left(\rho \right) \ \ [{\rm g \ cm}^{-3}]$', extend='both')

        plt.savefig(plot_dir + str(n) + '_FullFrameSnapshot_rho.png', bbox_inches='tight', dpi=300)
        #plt.show()
        plt.close()

        print('plot ' + str(n) + ' done')
        n += 1

    print('All output files from ' + str(start_nr) + ' to ' + str(n - 1) + ' have been plotted.' )

#######################################################################################################################

# define input and parameters
f = '2D-test'                                    # name of folder
dir = '/Users/wiebe/Documents/Running-Athena/Reading_data/' + f + '/'       # path to folder
plot_dir = dir + 'plots/'                                   # save plots here

a = 7.47e11        # orbital distance in cm -- a = 0.05 au
rp = 6.71e9        # planetary radius in cm -- WASP-107b

SCALE = 1          # global scaling parameter for pressure and radius

appendix = 'pws_3D'
orb = pw.read_trackfile(dir + 'pm_trackfile.dat')

start_nr = 0        # selection of snapshot files in the directory, if only one snapshot: 0, 1
end_nr = 3

mylevel = 2         # resolution level, choose between 1 and 5


plot_2D_frame(dir, plot_dir, start_nr, end_nr, xy=True)

