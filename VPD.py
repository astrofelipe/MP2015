from __future__ import division, print_function
from sklearn import mixture
from astroML.density_estimation import XDGMM
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse

#ARGUMENTOS
parser = argparse.ArgumentParser(description='VPD Plot')
parser.add_argument('<Input List>', help='Catalogo final de PMs')
parser.add_argument('--max-mag', type=float, default=20.0, help='Corte superior en magnitud (Default 20)')
parser.add_argument('--min-mag', type=float, default=8.0, help='Corte inferior en magnitud (Default 8)')
parser.add_argument('--max-err', type=float, default=1.0, help='Maximo error a considerar (Default 1)')
parser.add_argument('--comp', type=int, default=3, help='Nro de componentes para el Gaussian Mixture (Default 3)')
parser.add_argument('--center', nargs=2, default=None, help='Forzar centro a las coordenadas entregadas')
parser.add_argument('--hexbins', type=int, default=None, help='Usa bines hexagonales, se debe especificar tamano grilla')
parser.add_argument('--no-save', action='store_true', help='Mostrar plot en pantalla en vez de guardar')
parser.add_argument('--output', type=str, default='VPDbins.png', help='Cambiar nombre del file de output')

args = parser.parse_args()

inputs  = vars(args)['<Input List>']
max_mag = args.max_mag
min_mag = args.min_mag
max_err = args.max_err
n_comp  = args.comp

#Lee datos
#ID L B Ks EKs H EH J EJ Y EY Z EZ ULCOSB EULCOSB UB EUB NFRAMES
ids, magK, pmx, pmex, pmy, pmey, nframes = np.genfromtxt(inputs, unpack=True, usecols=(0,3,13,14,15,16, 17))

#Filtros
mag_mask = (magK < max_mag) & (magK > min_mag)
err_mask = (pmex**2 + pmey**2)**0.5 < max_err
nfr_mask = (nframes >= 8)
pm_mask  = (pmx**2 + pmy**2)**0.5 < 30
mask     = mag_mask & err_mask & nfr_mask & pm_mask

data     = np.transpose([pmx, pmy])[mask]
data_err = np.zeros(data.shape + data.shape[-1:])
diag     = np.arange(data.shape[-1])
data_err[:, diag, diag] = np.vstack([pmex**2, pmey**2]).T[mask]

print('\tTotal de estrellas:                %d' % len(mask))
print('\tNumero de estrellas seleccionadas: %d' % mask.sum())

#Calcula centros
g    = mixture.GMM(n_components=n_comp).fit(data)
x, y = np.transpose(g.means_)
print(g.means_)
idx  = np.argmin((x**2 + y**2)**0.5)

radio = ((pmx[mask] - x[idx])**2 + (pmy[mask] - y[idx])**2)**0.5 < 3
print(radio.sum())

#Plot
bins = np.arange(-15,15,1)
fig, ax = plt.subplots(figsize=[10,10])

if args.hexbins != None:
    h = ax.hexbin(data.T[0], data.T[1], gridsize=args.hexbins)
else:
    H, xedges, yedges, img = ax.hist2d(data.T[0], data.T[1], bins=bins)
    extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
    h = ax.matshow(np.rot90(H), extent=extent)
    ax.grid(linestyle='-', color='white', lw=.5, which='both', alpha=.2)
    ax.minorticks_on()

ax.plot(x, y, 'ow')
div = make_axes_locatable(ax)
cax = div.append_axes("right", size="5%", pad=0.2)

fig.colorbar(h, cax=cax)
ax.set_aspect('equal')

if args.no_save:
    plt.show()
else:
    fig.savefig(args.output, dpi=200, bbox_inches='tight')

import sys
sys.exit(1)
