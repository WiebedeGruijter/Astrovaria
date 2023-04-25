import numpy as np
import matplotlib.pyplot as plt
import imageio
import glob2


# define input and parameters
f = '0.05_DN_2D_V8_l5_100_e0.20'                                           # name of folder
dvar = 'div'
dir = '/Users/wiebe/Documents/Running-Athena/' + f + '/plots/' + dvar + '/'       # path to folder

list = []
for file_dir in glob2.glob(dir + '*.png'):
    list.append(file_dir)

list.sort()

start_nr = 50 # selection of snapshot files in the directory, if only one snapshot: 0, 1
end_nr= 151
list = list[start_nr:end_nr]
print(list)


# build gif
with imageio.get_writer(dir+'mygif_'+dvar+'.gif', mode='I') as writer:
    for filename in list:
        image = imageio.imread(filename)
        writer.append_data(image)
        print(str(filename) + ' appended')

print('gif has been plotted')
