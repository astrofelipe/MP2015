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
parser.add_argument('--lim', type=float, default=16, help='Limite en PM para el plot (cuadrado)')
parser.add_argument('--pm-cut', type=float, default=30, help='No considerar PM mayor al valor (Default 30)')
parser.add_argument('--comp', type=int, default=2, help='Nro de componentes para el Gaussian Mixture (Default 2)')
parser.add_argument('--center', nargs=2, default=None, help='Forzar centro a las coordenadas entregadas')
parser.add_argument('--hexbins', type=int, default=None, help='Usa bines hexagonales, se debe especificar tamano grilla')
parser.add_argument('--levels', type=int, default=7, help='Numero de niveles para el contour plot')
parser.add_argument('--hist2d', action='store_true', help='Hace el histograma en 2D en vez del KDE')
parser.add_argument('--no-save', action='store_true', help='Mostrar plot en pantalla en vez de guardar')
parser.add_argument('--output', type=str, default='VPDbins.png', help='Cambiar nombre del file de output')

args = parser.parse_args()

inputs  = vars(args)['<Input List>']
max_mag = args.max_mag
min_mag = args.min_mag
max_err = args.max_err
n_comp  = args.comp
lim     = args.lim

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
pm_mask  = (pmx**2 + pmy**2)**0.5 < args.pm_cut
mask     = mag_mask & err_mask & nfr_mask & pm_mask

data     = np.transpose([pmx, pmy])[mask]
data_err = np.zeros(data.shape + data.shape[-1:])
diag     = np.arange(data.shape[-1])
data_err[:, diag, diag] = np.vstack([pmex**2, pmey**2]).T[mask]

print('\tTotal de estrellas:                %d' % len(mask))
print('\tNumero de estrellas seleccionadas: %d' % mask.sum())
print('\tEstrellas con PM sobre %.2f: %d' % (args.pm_cut, (~pm_mask).sum()))

#Calcula centros
g    = mixture.GMM(n_components=n_comp, covariance_type='full').fit(data)
x, y = np.transpose(g.means_)
print('\nX_3D Y_3D    (cruces blancas)')
modulo = np.sum(g.means_**2, axis=1)**0.5
ormod  = np.argsort(modulo)
g.means_ = g.means_[ormod]
#g.means_[0] = np.array([0,0])
print(g.means_)

print('sig_X_3D sig_Y_3D')
for cov in g.covars_:
    print(np.sqrt(np.diag(cov)))

idx  = np.argmin((x**2 + y**2)**0.5)

##XDGMM
#XD = XDGMM(n_components=n_comp, n_iter=10).fit(data, data_err)
#xd_x, xd_y = np.transpose(XD.mu)
#print('X_XD Y_XD')
#print(XD.mu)

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
bins  = np.arange(-lim,lim+0.5,0.5)
cmap  = cm.get_cmap('jet')
color = cmap(np.linspace(0, 1, cmap.N))

fig = plt.figure(figsize=[10,10])

gs  = gridspec.GridSpec(4,4)

ax = plt.subplot(gs[1:,:-1])
axu = plt.subplot(gs[0,:-1])
axd = plt.subplot(gs[1:,-1])

#KDE
kde = KernelDensity(kernel='gaussian').fit((pmx[mask])[:,np.newaxis])
xx  = np.arange(-lim, lim, 0.05)
yy  = np.exp(kde.score_samples(xx[:,np.newaxis]))
x2d = xx[np.argmax(yy)]

a0  = np.exp(kde.score_samples(x[:,np.newaxis]))

if args.comp == 3:
    p0  = [a0[0]/2.0, x[0], 3, a0[1]/2.0, x[1], 3, a0[2]/2.0, x[2], 3]
elif args.comp == 2:
    p0  = [a0[0]/2.0, x[0], 3, a0[1]/2.0, x[1], 3]
elif args.comp == 1:
    p0  = [a0[0]/2.0, x[0], 3]

popt, pcov = curve_fit(gf, xx, yy, p0=p0, maxfev=100000)
x0g = popt[1::3]
x0s = popt[2::3]
x0e = np.sqrt(np.diag(pcov)[1::3])

axu.hist(pmx[mask], bins=bins, histtype='stepfilled', normed=True, color=cmap(0), alpha=.75)
axu.plot(xx, gaussian(xx, *popt[0:3]), color=cmap(0.6), lw=2.5)
if args.comp >= 2:
    axu.plot(xx, gaussian(xx, *popt[3:6]), color=cmap(0.4), lw=2.5)
if args.comp >= 3:
    axu.plot(xx, gaussian(xx, *popt[6:9]), color=cmap(0.2), lw=2)
axu.plot(xx,yy, color=cmap(0.8), lw=1.5, alpha=.75)

kde = KernelDensity(kernel='gaussian').fit((pmy[mask])[:,np.newaxis])
xx  = np.arange(-lim, lim, 0.05)
yy  = np.exp(kde.score_samples(xx[:,np.newaxis]))
y2d = xx[np.argmax(yy)]

a0  = np.exp(kde.score_samples(y[:,np.newaxis]))
if args.comp == 3:
    p0  = [a0[0]/2.0, y[0], 3, a0[1]/2.0, y[1], 3, a0[2]/2.0, y[2], 3]
elif args.comp == 2:
    p0  = [a0[0]/2.0, y[0], 3, a0[1]/2.0, y[1], 3]
elif args.comp == 1:
    p0  = [a0[0]/2.0, x[0], 3]

popt, pcov = curve_fit(gf, xx, yy, p0=p0, maxfev=100000)
y0g = popt[1::3]
y0s = popt[2::3]
y0e = np.sqrt(np.diag(pcov)[1::3])

axd.hist(pmy[mask], bins=bins, histtype='stepfilled', normed=True, color=cmap(0), alpha=.75, orientation='horizontal')
axd.plot(gaussian(xx, *popt[0:3]), xx, color=cmap(0.6), lw=2.5)
if args.comp >= 2:
    axd.plot(gaussian(xx, *popt[3:6]), xx, color=cmap(0.4), lw=2.5)
if args.comp >= 3:
    axd.plot(gaussian(xx, *popt[6:9]), xx, color=cmap(0.2), lw=2)
axd.plot(yy,xx, color=cmap(0.8), lw=1.5, alpha=.75)

print('\nX_2D Y_2D    (puntos blancos)')
print(np.transpose([x0g, y0g]))
print('X_err_2D Y_err_2D')
print(np.transpose([x0e, y0e]))
print('sig_X_2D sig_Y_2D')
print(np.transpose([x0s, y0s]))

print('\nMax 2D    (cuadrado blanco)')
print(x2d, y2d)

#print('\nCalculando KDE 2D...')
xy  = np.transpose([pmx, pmy])[mask]
k2d = KernelDensity(kernel='gaussian', bandwidth=0.4).fit(xy)

xg   = np.linspace(-lim, lim, 100)
X, Y = np.meshgrid(xg, xg)
XY   = np.vstack([Y.ravel(), X.ravel()]).T
logd = k2d.score_samples(XY)
Z    = np.exp(logd).reshape(X.shape)

x3d  = Y.ravel()[np.argmax(Z)]
y3d  = X.ravel()[np.argmax(Z)]
print('\nMax 3D    (triangulo blanco)')
print(x3d, y3d)

if args.hexbins != None:
    h = ax.hexbin(data.T[0], data.T[1], gridsize=args.hexbins)

elif args.hist2d:
    H, xedges, yedges, img = ax.hist2d(data.T[0], data.T[1], bins=bins)
    extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
    h = ax.matshow(np.rot90(H), cmap=cmap,extent=extent)
    ax.grid(linestyle='-', color='white', lw=.5, which='both', alpha=.2)
    ax.minorticks_on()

else:
    #h = ax.matshow(np.rot90(Z), extent=[-15, 15, -15, 15])
    h = ax.contourf(Y,X,Z)
    h = ax.contourf(Y,X,Z, levels=np.linspace(0, h.levels[-1], args.levels+1))
    ax.minorticks_on()

ax.plot([0], [0], '+', color='k', ms=15, mew=1.5)
ax.plot(x, y, 'xw', mew=1.5)
#ax.plot(xd_x, xd_y, 'x', color='gray', mew=1.5)
ax.plot(x0g, y0g, 'ow')
ax.plot(x2d, y2d, 'sw', ms=3)
ax.plot(x3d, y3d, '^w', ms=7)
ax.set_xlim(-lim, lim)
div = make_axes_locatable(ax)
#cax = div.append_axes("right", size="5%", pad=0.2)

co = fig.add_axes([0.29, 0.16, 0.25, 0.01])
cb = fig.colorbar(h, cax=co, orientation='horizontal')
#ax.add_patch(patches.Rectangle((-10, -22.1), 20, 1.5, alpha=.66, color='w', lw=0))
cb.ax.tick_params(labelcolor='#000000', pad=7.5, length=4)
plt.setp(plt.xticks()[1], rotation=45)
#ax.set_aspect('equal')

if args.no_save:
    plt.show()
else:
    fig.savefig(args.output, dpi=100, bbox_inches='tight')

import sys
sys.exit(1)
