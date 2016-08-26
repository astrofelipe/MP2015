import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import cm
from matplotlib.ticker import FormatStrFormatter
from matplotlib import patches
import argparse

#ARGUMENTOS
parser = argparse.ArgumentParser(description='VPD Plot')
parser.add_argument('Input', help='Catalogo final -> Fotometria + PMs')
parser.add_argument('--max-mag', type=float, default=20.0, help='Corte superior en magnitud (Default 20)')
parser.add_argument('--min-mag', type=float, default=5.0, help='Corte inferior en magnitud (Default 5)')
parser.add_argument('--max-err', type=float, default=10.0, help='Maximo error a considerar (Default 10)')
parser.add_argument('--lim', type=float, default=16, help='Limite en PM para el plot (cuadrado)')
parser.add_argument('--labelx', type=str, default='\mu_\ell \cos(b)', help='x axis label')
parser.add_argument('--labely', type=str, default='\mu_b', help='y axis label')
parser.add_argument('--invertx', action='store_true', help='Invierte el eje horizontal')
parser.add_argument('--rotation-right', type=int, default=0, help='Rotacion de los ticks del panel derecho')
parser.add_argument('--no-save', action='store_true', help='Mostrar plot en pantalla en vez de guardar')
parser.add_argument('--output', type=str, default='VPDfigure.png', help='Nombre del archivo de salida')

args = parser.parse_args()
bins = np.arange(-args.lim, args.lim + 0.5, 0.5)
lim  = args.lim
labx = args.labelx
laby = args.labely

#Archivo
#ID L B Ks EKs H EH J EJ Y EY Z EZ ULCOSB EULCOSB UB EUB NFRAMES
data = np.genfromtxt(args.Input, usecols=(3,13,14,15,16))

#Filtros
magf = (data[:,0] > args.min_mag) & (data[:, 0] < args.max_mag)
errf = data[:,2]**2 + data[:,4]**2 < args.max_err**2
mask = magf & errf

magK, pmx, pmex, pmy, pmey = data[mask].T

#Config plot
plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)
plt.rc('axes', labelsize=20)

#Plot
fig = plt.figure(figsize=[10,10])

hess = plt.axes([0.1, 0.1, 0.65, 0.65 ])
hx   = plt.axes([0.1, 0.1 + 0.65 + 0.02, 0.65, 0.2])
hy   = plt.axes([0.1 + 0.65 + 0.02, 0.1, 0.2, 0.65])

cmap = cm.get_cmap('RdYlBu_r')
cmap.set_bad('white', 1)

#Hess
H, xedges, yedges, img = hess.hist2d(pmx, pmy, bins=bins, cmap=cmap)
extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
#H      = np.ma.masked_where(H==0, H)
#im     = hess.matshow(np.flipud(np.rot90(H)), cmap=cmap, extent=extent, origin='lower')
im     = hess.pcolormesh(xedges, yedges, np.flipud(np.rot90(H)), cmap=cmap)#, extent=extent, origin='lower')
hess.xaxis.tick_bottom()
hess.xaxis.set_ticks_position('both')
hess.set_xlim(-lim, lim)
hess.set_ylim(-lim, lim)
hess.set_xlabel('$%s\quad\mathrm{(mas\ yr^{-1})}$' % labx)
hess.set_ylabel('$%s\quad\mathrm{(mas\ yr^{-1})}$' % laby)
hess.xaxis.set_major_formatter(FormatStrFormatter(r'$\mathrm{%d}$'))
hess.yaxis.set_major_formatter(FormatStrFormatter(r'$\mathrm{%d}$'))

#ax.plot(x,y, alpha=.2)

#Colorbar
co = fig.add_axes([0.25, 0.16, 0.35, 0.015])
cb = fig.colorbar(im, cax=co, orientation='horizontal', format=FormatStrFormatter(r'$\mathrm{%d}$'))
cb.ax.tick_params(labelcolor='w')

#Histograma PMX
hx.hist(pmx, bins=bins, histtype='stepfilled', color=cmap(0.05))
hx.xaxis.set_ticks_position('top')
hx.set_yticks(hx.get_yticks()[::2])
hx.set_xlim(-lim, lim)
hx.set_ylabel('$N$')
hx.xaxis.set_major_formatter(FormatStrFormatter(r'$\mathrm{%d}$'))
hx.yaxis.set_major_formatter(FormatStrFormatter(r'$\mathrm{%d}$'))

#Histogama PMY
hy.hist(pmy, bins=bins, histtype='stepfilled', color=cmap(0.05), orientation='horizontal')
hy.set_ylim(-lim, lim)
hy.yaxis.set_ticks_position('right')
hy.set_xticks(hy.get_xticks()[::2])
plt.setp(hy.get_xticklabels(), rotation=args.rotation_right)
hy.set_xlabel('$N$')
hy.xaxis.set_major_formatter(FormatStrFormatter(r'$\mathrm{%d}$'))
hy.yaxis.set_major_formatter(FormatStrFormatter(r'$\mathrm{%d}$'))

hess.set_aspect('equal')
hess.tick_params(which='both', length=5, color='w')

if args.invertx:
    hx.invert_xaxis()
    hess.invert_xaxis()

if args.no_save:
    plt.show()
else:
    fig.savefig(args.output, dpi=200, bbox_inches='tight')
