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
        orb = pw.read_trackfile(dir+filename+'/' + 'pm_trackfile.dat')
   
        list = []
        for file_dir in glob2.glob(dir+filename+'/' + '*.athdf'):
            list.append(file_dir)
        list.sort()

        list = list[start_nr:end_nr]
        tot_list.append(list)

    print(tot_list)

    n = 0 + start_nr
    for i in range(len(list)):

        fig, ax = plt.subplots(1, len(tot_list), figsize=((15,15)), sharey=True)
        lim = 2e12
        filename = dvar
        
        for j in range(len(tot_list)):
            snapshot = tot_list[j][i]

            # read snapshot
            dblank = pw.ar.athdf(snapshot, quantities=[], level=mylevel, subsample=True)

            if xy == True:
                #x1_max = r[cm], x2_max = theta[rad], x3_max = phi[rad]
                x2sliceval = dblank['x2v'][np.argmin(np.abs(dblank['x2v'] - np.pi / 2))]
                d = pw.read_data(snapshot, orb, level=mylevel, x2_max=x2sliceval, x2_min=x2sliceval)
                X = pw.get_plot_array_midplane(d['x'][:, 0, :]) / lenscale
                Y = pw.get_plot_array_midplane(d['y'][:, 0, :]) / lenscale
                Z = pw.get_plot_array_midplane(np.log10(d[dvar][:, 0, :] * SCALE))
                filename = dvar+'_xy'
            else:
                # gives only left side of the full frame in xz plane
                x3sliceval = dblank['x3v'][np.argmin(np.abs(dblank['x3v']- np.pi))]
                d = pw.read_data(snapshot, orb, level=mylevel, x3_min=x3sliceval, x3_max=x3sliceval)
                X = pw.get_plot_array_midplane(d['x'][0, :, :]) / lenscale
                Y = pw.get_plot_array_midplane(d['z'][0, :, :]) / lenscale
                Z = pw.get_plot_array_midplane(np.log10(d[dvar][0, :, :] * SCALE))
                filename = dvar+'_xz'

            
            if zoom_level!=None:
                im = ax[j].pcolormesh(X+a/lenscale, Y, Z, cmap='magma', vmin=vmin, vmax=vmax, shading='nearest', rasterized=True)
                lim = (zoom_level*lenscale)
                filename = filename+'_zoomed'
            else:
                im = ax[j].pcolormesh(X, Y, Z, cmap='magma', vmin=vmin, vmax=vmax, rasterized=True)
            
            ax[j].set_aspect('equal')
            ax[j].set_xlim(-lim / lenscale, lim / lenscale)
            ax[j].set_ylim(-lim / lenscale, lim / lenscale)
            ax[j].set_title('e = 0.'+f[j][-2:])

            # labels and colorbar

            ax[j].set_xlabel(r'$x$ '+unit_name)
            if xy == True:
                ax[0].set_ylabel(r'$y$ '+unit_name)
            else: 
                ax[0].set_ylabel(r'$z$ '+unit_name)

        cax = fig.add_axes([0.96, 0.43, 0.01, 0.15])
        if dvar == 'rho':
            cb = fig.colorbar(im, cax=cax, label=r'$\log_{10}\left(\rho \right) \ \ [{\rm g \ cm}^{-3}]$', extend='both')
            
        else:
            cb = fig.colorbar(im, cax=cax, label=r'$\log_{10}\left({\rm P} \right) \ \ [{\rm barye}]$', extend='both')

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

plot_dir = dir + 'multiplots/'                                   # save plots here

a = 7.47e11      # orbital distance in cm -- a = 0.05 au
rp = 6.71e9        # planetary radius in cm -- WASP-107b

SCALE = 1          # global scaling parameter for pressure and radius

appendix = 'pws_3D'

start_nr = 0        # selection of snapshot files in the directory, if only one snapshot: 0, 1
end_nr = 301

mylevel = 3         # resolution level, choose between 1 and 5

#Variables:     density: 'rho' colormap: 'magma'
#               pressure: 'press' colormap: 'plasma'


plot_2D_frame(dir, plot_dir, start_nr, end_nr, xy=True, r_p=True, dvar='rho', zoom_level=None, vmin=None, vmax=None)

