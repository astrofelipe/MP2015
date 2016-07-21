import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
import argparse
from matplotlib import rc
import numpy as np
import sys

#README
parser = argparse.ArgumentParser(description='CMD Bins')
parser.add_argument('INPUT', help='Catalogo final de PMs')
parser.add_argument('x0', type=float, help='Coordenada X para elegir estrellas')
parser.add_argument('y0', type=float, help='Coordenada Y para elegir estrellas')
parser.add_argument('--filter', type=str, default='J', help='Segundo filtro para confeccionar color (Default J)')
parser.add_argument('--mags', type=float, nargs=2, default=(10, 18), help='Limites en magnitud (Default 10 18)', metavar=('min', 'max'))
parser.add_argument('--dmag', type=float, default=1.0, help='Intervalo para hacer los cortes (Default 1)', metavar='dmag')
parser.add_argument('--max-err', type=float, default=10.0, help='Maximo error a considerar (Default 10)', metavar='err')
parser.add_argument('--color', type=float, nargs=2, default=(0, 2), help='Limite en color', metavar=('min', 'max'))
parser.add_argument('--radios', type=float, nargs=2, default=(4, 9.5), help='Limite en radios para seleccionar (Default 4 9.5)', metavar=('min', 'max'))
parser.add_argument('--dradio', type=float, default=0.5, help='Incremento en radio de seleccion por bin de magnitud (Default 0.5)', metavar='dr')
parser.add_argument('--lim', type=float, nargs=2, default=(-24, 24), help='Limite en los ejes del VPD', metavar=('min', 'max'))
parser.add_argument('--no-save', action='store_true', help='Mostrar plot en pantalla en vez de guardar')

args = parser.parse_args()

##############
# PARAMETROS #
##############
file1 = args.INPUT

b2nam            = args.filter
min_mag, max_mag = args.mags
dmag             = args.dmag
col_min, col_max = args.color
max_err          = args.max_err
rad_pm           = np.arange(args.radios[0], args.radios[1], args.dradio)
pm_min, pm_max   = args.lim


x0, y0  = args.x0, args.y0
#x0, y0  = 10.52710485, 6.82687392  #Centro para la seleccion

#Plots
rc('text', usetex=True) #Usa TeX para las fuentes
rc('xtick', labelsize=15) #Tamano fuente ticks en x
rc('ytick', labelsize=15) #Tamano fuente ticks en y
rc('xtick.major', size=7.5) #Tamano fuente ticks en x
rc('ytick.major', size=7.5) #Tamano fuente ticks en x
rc('axes', labelsize=24) #Tamano fuente labels

cmd_psize = 3   #Tamano puntos CMD
pmd_psize = 2   #Tamano puntos diagrama PM

pm_ticks  = np.arange(-30,31,15) #Ticks para ejes del diagrama de PM. Intervalo de 15 en el rango [-30,31) para que incluya al 30
mag_ticks = np.arange(10,20,1)   #Ticks para eje de magnitud. Intervalo de 1, desde 9 a 20. Son solo los ticks, los limites mandan en el resultado final.

#Esto hace un circulo nada mas
def circle(x0, y0, r):
    theta = np.linspace(0, 2*np.pi, 1000)
    x     = x0 + r*np.cos(theta)
    y     = y0 + r*np.sin(theta)
    return x,y


#########################
# LEE ARCHIVOS Y FILTRA #
#########################
#PM_final
data = np.genfromtxt(file1, usecols=(0,3,7,13,14,15,16))
#Filtra por error
pme = (data[:,4]**2 + data[:,6]**2)**0.5
maskerr = pme <= max_err

ids, magK, magB, pmx, pmex, pmy, pmey = data[maskerr].T
ids = ids.astype(int)

########
# PLOT #
########
mags = np.arange(min_mag, max_mag+1, dmag)

#Tweak
#mags[-1] = 18.2

nint = len(mags) - 1
bands = ['K_s', 'H', 'J', 'Y' ,'Z']

fig  = plt.figure(figsize=[5*5.5/2, nint*3/2])
grid = gs.GridSpec(nint, 5)

pm = [plt.subplot(grid[i,0]) for i in range(nint)]
to = plt.subplot(grid[:,1:3])
fi = plt.subplot(grid[:,3:])

filtro = np.zeros(len(pmx)).astype(bool)

#Itera haciendo los cortes en magnitud para el plot de PMs
for i in range(nint):
    maskmag = (magK < mags[i+1]) & (magK > mags[i])

    x, y = pmx[maskmag], pmy[maskmag]
    pm[i].plot(x, y, '.k', ms=pmd_psize, rasterized=True)

    #Seleccion
    maskred = ((pmx[maskmag] - x0)**2 + (pmy[maskmag] - y0)**2)**0.5 <= rad_pm[i]
    filtro[maskmag] += maskred
    xc, yc  = circle(x0, y0, rad_pm[i])

    xr, yr  = x[maskred], y[maskred]
    pm[i].plot(xr, yr, '.r', ms=pmd_psize)#, rasterized=True)
    pm[i].plot(xc, yc, '-', c='lime', lw=1.5)

    pm[i].set(xlim=(pm_min, pm_max), ylim=(pm_min, pm_max), yticks=pm_ticks, xticks=pm_ticks, aspect='equal', adjustable='box')

    if i == int(nint/2):
        if nint%2 == 0:
            pm[i].set_ylabel('$\_ \quad\quad\quad\ \mu_b$  $\mathrm{mas\ yr^{-1}}$')
        else:
            pm[i].set_ylabel('$\mu_b$  $\mathrm{mas\ yr^{-1}}$')


    #Plotea en el panel derecho
    to.plot((magB-magK)[maskmag][maskred], magK[maskmag][maskred], '.k', ms=cmd_psize)#, rasterized=True)
    fi.plot((magB-magK)[maskmag][~maskred], magK[maskmag][~maskred], '.k', ms=cmd_psize)#, rasterized=True)

pm[-1].set_xlabel('$\mu_\ell\cos b$  $\mathrm{mas\ yr^{-1}}$')

#to.plot(magB-magK, magK, '.k', ms=cmd_psize, rasterized=True)
#min_mag -= 0.2
#max_mag += 0.2  #Tweak

to.set(xlim=(col_min, col_max), ylim=(max_mag, min_mag), ylabel='$K_s$', xlabel='$%s - K_s$' % b2nam, yticks=mag_ticks)
fi.set(xlim=(col_min, col_max), ylim=(max_mag, min_mag), ylabel='$K_s$', xlabel='$%s - K_s$' % b2nam, yticks=mag_ticks)

#Guarda catalogos
fmt  = '%d %.6f %.6f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.3f %.2f %.3f %.2f %d'
hdr  = 'ID L B Ks EKs H EH J EJ Y EY Z EZ ULCOSB EULCOSB UB EUB NFRAMES'
data = np.genfromtxt(file1)

data_in  = data[maskerr][filtro]
data_out = data[maskerr][~filtro]
np.savetxt(file1.replace('.gale', '_in.gale'), data_in, fmt=fmt, header=hdr)
np.savetxt(file1.replace('.gale', '_out.gale'), data_out, fmt=fmt, header=hdr)


#Cambios finales
grid.update(hspace=0, wspace=0.8)            #Espacio vertical nulo, horizontal mayor
if args.no_save:
    plt.show()
else:
    fig.savefig(file1.replace('.gale', '.png'), dpi=200, bbox_inches='tight')
