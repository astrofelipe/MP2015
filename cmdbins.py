import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
from matplotlib import rc
import numpy as np
import sys
##############
# PARAMETROS #
##############
#Segundo filtro
band2 = 7 #Para confeccionar color (band2 - K) -- K, H, J, Y, Z vienen en las columnas 3, 5, 7, 9, 11
file1 = sys.argv[1]    #Archivo con los PM
file2 = sys.argv[2]    #Archivo para extraer el filtro

#Magnitudes
min_mag = 10.0
max_mag = 18.0  #Limites en magnitud para el plot
dmag    = 1.0   #Intervalo para hacer los cortes

col_min = 0.0
col_max = 2.0     #Limites en color (eje x)

#PMs
max_err = 10.0   #Error maximo (modulo) para puntos a mostrar
pm_magl = 15    #Estrellas con magnitud mayor usan rad_pm2
rad_pm  = np.arange(5,9.5,0.5) #Radios para seleccionar estrellas (arreglo largo nro de bines)

x0, y0  = -6.82687392, 10.52710485  #Centro para la seleccion

pm_min = -24
pm_max = 24     #Limites en PMX y PMY

#Plots
rc('text', usetex=True) #Usa TeX para las fuentes
rc('xtick', labelsize=15) #Tamano fuente ticks en x
rc('ytick', labelsize=15) #Tamano fuente ticks en y
rc('xtick.major', size=7.5) #Tamano fuente ticks en x
rc('ytick.major', size=7.5) #Tamano fuente ticks en x
rc('axes', labelsize=24) #Tamano fuente labels

cmd_psize = 3   #Tamano puntos CMD
pmd_psize = 1   #Tamano puntos diagrama PM

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
data = np.genfromtxt(file1, usecols=(0,3,4,5,8,9))
#Filtra por error
pme = (data[:,4]**2 + data[:,5]**2)**0.5
maskerr = pme <= max_err

ids, pmx, pmy, magK, pmex, pmey = data[maskerr].T
ids = ids.astype(int)

#Segundo archivo
id2, mag2 = np.genfromtxt(file2, unpack=True, usecols=(0,band2), skip_header=3)
id2 = id2.astype(int)

orden = np.argsort(id2)
id2   = id2[orden]
mag2  = mag2[orden]
mag2[mag2>98] = np.nan

#Descarta los IDs 5000000 (sin fotometria en Ks)
#mask5 = id2 < 5000000
#id2   = id2[mask5]
#mag2  = mag2[mask5]

#Interseccion y magnitud en banda 2 a usar
magB    = np.empty_like(magK)
magB[:] = np.nan

inter1 = np.in1d(ids, id2)
inter2 = np.in1d(id2, ids)
magB[inter1]  = mag2[inter2]

########
# PLOT #
########
mags = np.arange(min_mag, max_mag+1, dmag)

#Tweak
mags[-1] = 18.2

nint = len(mags) - 1
bands = ['K_s', 'H', 'J', 'Y' ,'Z']
b2nam = bands[int((band2-1)/2 - 1)]

fig  = plt.figure(figsize=[5*5.5/2, nint*3/2])
grid = gs.GridSpec(nint, 5)

pm = [plt.subplot(grid[i,0]) for i in range(nint)]
to = plt.subplot(grid[:,1:3])
fi = plt.subplot(grid[:,3:])

#Itera haciendo los cortes en magnitud para el plot de PMs
for i in range(nint):
    maskmag = (magK < mags[i+1]) & (magK > mags[i])

    x, y = pmx[maskmag], pmy[maskmag]
    pm[i].plot(x, y, '.k', ms=pmd_psize, rasterized=True)

    #Seleccion
    maskred = ((pmx[maskmag] - x0)**2 + (pmy[maskmag] - y0)**2)**0.5 <= rad_pm[i]
    xc, yc  = circle(x0, y0, rad_pm[i])
    #if mags[i+1]-dmag < pm_magl:
    #    maskred = ((pmx[maskmag] - x0)**2 + (pmy[maskmag] - y0)**2)**0.5 <= rad_pm1
    #    xc, yc  = circle(x0, y0, rad_pm1)
    #else:
    #    maskred = ((pmx[maskmag] - x0)**2 + (pmy[maskmag] - y0)**2)**0.5 <= rad_pm2
    #    xc, yc  = circle(x0, y0, rad_pm2)

    xr, yr  = x[maskred], y[maskred]
    pm[i].plot(xr, yr, '.r', ms=pmd_psize)#, rasterized=True)
    pm[i].plot(xc, yc, '-', c='lime', lw=1.5)

    pm[i].set(xlim=(pm_min, pm_max), ylim=(pm_min, pm_max), yticks=pm_ticks, xticks=pm_ticks, aspect='equal', adjustable='box')

    if i == int(nint/2):
        if nint%2 == 0:
            pm[i].set_ylabel('$\_ \quad\quad\quad\ \mu_y\enskip\mathrm{mas\ yr^{-1}}$')
        else:
            pm[i].set_ylabel('$\mu_y\enskip\mathrm{mas\ yr^{-1}}$')


    #Plotea en el panel derecho
    to.plot((magB-magK)[maskmag][maskred], magK[maskmag][maskred], '.k', ms=cmd_psize, rasterized=True)
    fi.plot((magB-magK)[maskmag][~maskred], magK[maskmag][~maskred], '.k', ms=cmd_psize, rasterized=True)

pm[-1].set_xlabel('$\mu_x\enskip\mathrm{mas\ yr^{-1}}$')

to.plot(magB-magK, magK, '.k', ms=cmd_psize, rasterized=True)
#min_mag -= 0.2
max_mag += 0.2  #Tweak

to.set(xlim=(col_min, col_max), ylim=(max_mag, min_mag), ylabel='$K_s$', xlabel='$%s - K_s$' % b2nam, yticks=mag_ticks)
fi.set(xlim=(col_min, col_max), ylim=(max_mag, min_mag), ylabel='$K_s$', xlabel='$%s - K_s$' % b2nam, yticks=mag_ticks)


#Cambios finales
grid.update(hspace=0, wspace=0.8)            #Espacio vertical nulo, horizontal mayor
fig.savefig('cmdcut.png', dpi=200, bbox_inches='tight')
