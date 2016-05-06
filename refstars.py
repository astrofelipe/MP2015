from __future__ import division, print_function

import matplotlib.pyplot as plt
import numpy as np
import sys
from matplotlib import patches
from matplotlib.path import Path
from astropy.utils.console import color_print

color_print('REFSTARS.py', 'yellow')

if len(sys.argv) == 1:
    color_print()

#Argumentos
cmd_file = sys.argv[1]  #CMD (.mat o .matgc)
ver_file = sys.argv[2]  #Vertices del poligono

data  = np.genfromtxt(cmd_file)
K     = data[:,5]
J_K   = data[:,12] - K
verts = np.genfromtxt(ver_file)

path = Path(verts)
idx  = path.contains_points(np.transpose([J_K,K]))

refstars = data[idx]
hdr = 'ID RA_1 DEC_1 X_1 Y_1 MAG_1 MAG_ERR_1 ID_2 RA_2 DEC_2 X_2 Y_2 MAG_2 MAG_ERR_2'
fmt = '%d %.7f %.7f %.3f %.3f %.3f %.3f %d %.7f %.7f %.3f %.3f %.3f %.3f'
np.savetxt('arefstars.dat', refstars, header=hdr, fmt=fmt)

'''
fig, ax = plt.subplots()
ax.plot(J_K, K, '.k', zorder=-100, ms=1.9)
ax.plot(J_K[idx], K[idx], '.g', ms=2)
#ax.plot(verts[:,0], verts[:,1], '.r')
ax.set_ylim(20,8)
plt.show()
'''
