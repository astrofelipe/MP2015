import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
import numpy as np
import matplotlib.patches as patches
import sys
from matplotlib import rc, cm
from matplotlib.colors import LinearSegmentedColormap

from astroML.density_estimation import XDGMM

#PARAMETROS
file1   = sys.argv[1]   #Archivo con los PM
pm_min  = -24
pm_max  = 24    #PM maximo y minimo a considerar
dbin    = 0.4   #Tamaño del bin

mag_min = 8
mag_max = 17    #Corte en magnitud

max_err = 5    #Error maximo a considerar en cada medicion

#Plots
rc('text', usetex=True) #Usa TeX para las fuentes
rc('xtick', labelsize=15) #Tamaño fuente ticks en x
rc('ytick', labelsize=15) #Tamaño fuente ticks en y
rc('xtick.major', size=7.5) #Tamaño fuente ticks en x
rc('ytick.major', size=7.5) #Tamaño fuente ticks en x
rc('axes', labelsize=24) #Tamaño fuente labels

pm_ticks  = np.arange(-30,31,10) #Ticks para ejes del diagrama de PM. Intervalo de 10 en el rango [-30,31) para que incluya al 30

#GAUSSIANAS
def gaussian(x, mu, sig):
    arg = -(x - mu)**2 / (2*sig**2)
    return np.exp(arg) / (sig*np.sqrt(2*np.pi))
#Suma de dos
def gauss2(x1, x2, mu1, mu2, sig1, sig2):
    return gaussian(x1, mu1, sig1) + gaussian(x2, mu2, sig2)
#Suma de tres
def gauss2(x1, x2, x3, mu1, mu2, mu3, sig1, sig2, sig3):
    return gaussian(x1, mu1, sig1) + gaussian(x2, mu2, sig2) + gaussian(x3, mu3, sig3)

#Carga los datos
data = np.genfromtxt(file1, usecols=(3,4,5,8,9))
#Filtra por error
pme     = (data[:,3]**2 + data[:,4]**2)**0.5
maskerr = pme <= max_err
maskmag = (data[:,2] > mag_min) & (data[:,2] < mag_max)

pmx, pmy, magK, pmex, pmey = data[maskerr*maskmag].T

########
# PLOT #
########
bins = np.arange(pm_min, pm_max + dbin, dbin)
cmap = cm.get_cmap('RdYlBu_r')
cmap2 = cm.get_cmap('viridis')
color = cmap2(np.linspace(0.33, 1, cmap.N // 2))
cmap  = LinearSegmentedColormap.from_list('Eeeee', color)

hist = plt.hist2d(pmx, pmy, bins=bins)
exte = [hist[2][0], hist[2][-1], hist[1][0], hist[1][-1]]

#Orienta el histograma de la forma correcta ([filas, columnas] = [y,x] no [x,y])
H = hist[0]
H = np.flipud(np.rot90(H))
H = np.ma.masked_where(H==0, H)

fig  = plt.figure(figsize=[4*3, 4*3])
grid = gs.GridSpec(4,4)

up = plt.subplot(grid[0,:3])
ri = plt.subplot(grid[1:,-1])
pm = plt.subplot(grid[1:,0:3])
#co = plt.subplot(grid[0,-1])

im = pm.imshow(H, cmap=cmap, extent=exte, origin='lower', interpolation='none')
pm.set(xlabel='$\mu_x\enskip\mathrm{mas\ yr^{-1}}$', ylabel='$\mu_y\enskip\mathrm{mas\ yr^{-1}}$', xlim=(pm_min, pm_max), ylim=(pm_min, pm_max), yticks=pm_ticks, xticks=pm_ticks, aspect='equal', adjustable='box')
pm.xaxis.set_ticks_position('bottom')
pm.tick_params(axis='x',which='major',top='on')

up.hist(pmx, histtype='stepfilled', bins=bins, lw=2, alpha=.9, facecolor=cmap(0.1), edgecolor=cmap(0.1))
up.set_xlim(pm_min, pm_max)
aspect = np.diff(up.get_xlim()) / np.diff(up.get_ylim())
up.set(ylabel='$N$', aspect=aspect[0]/3, adjustable='box')
up.xaxis.tick_top()
up.tick_params(axis='x',which='major',bottom='on')

ri.hist(pmy, histtype='stepfilled', bins=bins, lw=2, alpha=.9, facecolor=cmap(0.1), edgecolor=cmap(0.1), orientation='horizontal')
ri.set_ylim(pm_min, pm_max)
aspect = np.diff(ri.get_xlim()) / np.diff(ri.get_ylim())
ri.set(xlabel='$N$', aspect=aspect[0]*3.02, adjustable='box')
ri.yaxis.tick_right()
ri.tick_params(axis='y',which='major',left='on')
plt.setp(ri.get_xticklabels(), rotation=45)
#co = fig.add_axes([0.775, 0.70625, 0.05, 0.175])
co = fig.add_axes([0.29, 0.16, 0.25, 0.01])
#co.set_aspect(0.5)
#co.yaxis.tick_right()
#co.tick_params(axis='x', which='both', bottom='off')
cb = fig.colorbar(im, cax=co, orientation='horizontal')
pm.add_patch(patches.Rectangle((-10, -22.1), 20, 1.5, alpha=.66, color='w', lw=0))
cb.ax.tick_params(labelcolor='#000000', pad=7.5, length=4)
#cb.set_label('$N_{\mathrm{bin}}$')

grid.update(hspace=-0.1, wspace=0.01)
fig.savefig('hess.png', dpi=300, bbox_inches='tight')

#XDGMM
'''
print('XDGMM...')
dataset = np.vstack([pmx[magK<15], pmy[magK<15]]).T
dataerr = np.zeros(dataset.shape + dataset.shape[-1:])
diag    = np.arange(dataset.shape[-1])
dataerr[:, diag, diag] = np.vstack([pmex[magK<15]**2, pmey[magK<15]**2]).T

clf = XDGMM(3, 20, verbose=True)
clf.fit(dataset, dataerr)
print(clf.mu, clf.V, clf.alpha)
'''
