from __future__ import division, print_function
from sklearn import mixture
from sklearn.neighbors import KernelDensity
from astroML.density_estimation import XDGMM
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit
from astroML.plotting import hist
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import argparse

#ARGUMENTOS
parser = argparse.ArgumentParser(description='VPD Plot')
parser.add_argument('<Input List>', help='Catalogo final -> Fotometria + PMs')
parser.add_argument('--max-mag', type=float, default=20.0, help='Corte superior en magnitud (Default 20)')
parser.add_argument('--min-mag', type=float, default=8.0, help='Corte inferior en magnitud (Default 8)')
parser.add_argument('--max-err', type=float, default=2.0, help='Maximo error a considerar (Default 2)')
parser.add_argument('--comp', type=int, default=2, help='Nro de componentes para el Gaussian Mixture (Default 2)')
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
if 'PM_final' in inputs:
    ids, pmx, pmy, magK, nframes, pmex, pmey = np.genfromtxt(inputs, unpack=True, usecols=(0,3,4,5,6,8,9))
else:
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
g    = mixture.GMM(n_components=n_comp, covariance_type='full').fit(data)
x, y = np.transpose(g.means_)
print('X_guess Y_guess')
print(g.means_)
idx  = np.argmin((x**2 + y**2)**0.5)

#radio = ((pmx[mask] - x[idx])**2 + (pmy[mask] - y[idx])**2)**0.5 < 3

#Gaussianas
def gaussian(x, amp, mu, sig):
    if (amp < 0) or (sig < 0):
        return np.inf
    return amp * np.exp(-(x-mu)**2 / (2*sig**2))

def two_gaussian(x, amp1, mu1, sig1,
                    amp2, mu2, sig2):
    return (gaussian(x, amp1, mu1, sig1) +
            gaussian(x, amp2, mu2, sig2))

def three_gaussian(x, amp1, mu1, sig1,
                      amp2, mu2, sig2,
                      amp3, mu3, sig3):
    return (gaussian(x, amp1, mu1, sig1) +
    gaussian(x, amp2, mu2, sig2) +
    gaussian(x, amp3, mu3, sig3))

if args.comp == 2:
    gf = two_gaussian
elif args.comp == 3:
    gf = three_gaussian
elif args.comp == 1:
    gf = gaussian

#Plot
bins  = np.arange(-15,15,0.5)
cmap  = cm.get_cmap('jet')
color = cmap(np.linspace(0, 1, cmap.N))

fig = plt.figure(figsize=[10,10])

gs  = gridspec.GridSpec(4,4)

ax = plt.subplot(gs[1:,:-1])
axu = plt.subplot(gs[0,:-1])
axd = plt.subplot(gs[1:,-1])

#KDE
kde = KernelDensity(kernel='gaussian').fit((pmx[mask])[:,np.newaxis])
xx  = np.arange(-15, 15, 0.05)
yy  = np.exp(kde.score_samples(xx[:,np.newaxis]))

a0  = np.exp(kde.score_samples(x[:,np.newaxis]))

if args.comp == 3:
    p0  = [a0[0]/2.0, x[0], 3, a0[1]/2.0, x[1], 3, a0[2]/2.0, x[2], 3]
elif args.comp == 2:
    p0  = [a0[0]/2.0, x[0], 3, a0[1]/2.0, x[1], 3]
elif args.comp == 1:
    p0  = [a0[0]/2.0, x[0], 3]

popt, pcov = curve_fit(gf, xx, yy, p0=p0, maxfev=10000)
x0g = popt[1::3]
x0e = np.sqrt(np.diag(pcov)[1::3])

axu.hist(pmx[mask], bins=bins, histtype='stepfilled', normed=True, color=cmap(0), alpha=.75)
axu.plot(xx, gaussian(xx, *popt[0:3]), color=cmap(0.6), lw=2.5)
axu.plot(xx, gaussian(xx, *popt[3:6]), color=cmap(0.4), lw=2.5)
#axu.plot(xx, gaussian(xx, *popt[6:9]), color=cmap(0.2), lw=2)
axu.plot(xx,yy, color=cmap(0.8), lw=1.5, alpha=.75)

kde = KernelDensity(kernel='gaussian').fit((pmy[mask])[:,np.newaxis])
xx  = np.arange(-15, 15, 0.05)
yy  = np.exp(kde.score_samples(xx[:,np.newaxis]))

a0  = np.exp(kde.score_samples(y[:,np.newaxis]))
if args.comp == 3:
    p0  = [a0[0]/2.0, y[0], 3, a0[1]/2.0, y[1], 3, a0[2]/2.0, y[2], 3]
elif args.comp == 2:
    p0  = [a0[0]/2.0, y[0], 3, a0[1]/2.0, y[1], 3]
elif args.comp == 1:
    p0  = [a0[0]/2.0, x[0], 3]
    
popt, pcov = curve_fit(gf, xx, yy, p0=p0, maxfev=10000)
y0g = popt[1::3]
y0e = np.sqrt(np.diag(pcov)[1::3])

axd.hist(pmy[mask], bins=bins, histtype='stepfilled', normed=True, color=cmap(0), alpha=.75, orientation='horizontal')
axd.plot(gaussian(xx, *popt[0:3]), xx, color=cmap(0.6), lw=2.5)
axd.plot(gaussian(xx, *popt[3:6]), xx, color=cmap(0.4), lw=2.5)
#axd.plot(gaussian(xx, *popt[6:9]), xx, color=cmap(0.2), lw=2)
axd.plot(yy,xx, color=cmap(0.8), lw=1.5, alpha=.75)

print('X Y')
print(np.transpose([x0g, y0g]))
print('X_err Y_err')
print(np.transpose([x0e, y0e]))

if args.hexbins != None:
    h = ax.hexbin(data.T[0], data.T[1], gridsize=args.hexbins)
else:
    H, xedges, yedges, img = ax.hist2d(data.T[0], data.T[1], bins=bins)
    extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
    h = ax.matshow(np.rot90(H), extent=extent, cmap=cmap)
    ax.grid(linestyle='-', color='white', lw=.5, which='both', alpha=.2)
    ax.minorticks_on()

#ax.plot(x, y, 'og')
ax.plot(x0g, y0g, 'ow')
ax.plot([0], [0], '+', color='k', ms=15, mew=1.5)
div = make_axes_locatable(ax)
#cax = div.append_axes("right", size="5%", pad=0.2)

#fig.colorbar(h, cax=cax)
ax.set_aspect('equal')

if args.no_save:
    plt.show()
else:
    fig.savefig(args.output, dpi=200, bbox_inches='tight')

import sys
sys.exit(1)
