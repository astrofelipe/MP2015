import os
import sys
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn import linear_model
import multiprocessing
#from joblib import Parallel, delayed
from matplotlib import gridspec
from scipy.optimize import curve_fit
from astropy.io import ascii
from astropy.utils.console import ProgressBar
from astropy.table import Table, join, hstack
from pm_funcs import barra, linear_regression
from bayesian_regression import BayesianRegression
import pm_funcs

#PARAMETROS
nframes, nbins, limplotpm, nprocs = pm_funcs.get_pm1a1()
limplot = limplotpm

print 'Usando...'
print 'nframes: %d' % nframes
print 'nbins:   %d' % nbins

stilts_folder = os.path.dirname(os.path.realpath(__file__))
cpun = multiprocessing.cpu_count()

def recta(x,a,b):
    return a*x + b

def load_file(archivo):
    data = np.genfromtxt(archivo, unpack=True)
    return data

def load_table(archivo):
    tabla = Table.read(archivo, format='ascii')
    return tabla

#Lee los archivos con DX DY
referencia = sys.argv[1]

archivos   = np.sort(glob.glob('./PMs/PM_*'))
nro_arch   = len(archivos)
#threshold  = int(nro_arch * thr_per / 100)

nro_epoca = np.sort([int(f.split('_')[1].split('.')[0]) for f in archivos])
print 'Epocas: ', nro_epoca

#Realiza el match entre los PM_*.dat
if not os.path.isfile('PM.hdf5'):
    todos   = barra(load_file, archivos, nprocs)
    maximos = np.zeros(len(todos))
    refdatax = load_file(referencia)
    refdatax = refdatax.T[np.argsort(refdatax[0])].T

    for i in xrange(len(todos)):
        maximos[i] = np.max(todos[i][0])
    maximo = np.max(maximos)
    maximo = np.max([maximo, np.nanmax(refdatax[0])])

    total_id  = np.arange(1, maximo+1)

    refdata  = np.zeros((maximo, len(refdatax)))
    refdata[:] = np.nan


    refidsm = np.in1d(total_id, refdatax[0])
    refdata[refidsm] = refdatax.T

    #Genera matriz donde van todos los catalogos
    allcat    = np.zeros((maximo, len(todos)*4))
    allcat[:] = np.nan


    #Rellena la matriz
    print 'Ingresando datos a la matriz...'
    for i in ProgressBar(xrange(len(todos))):
        ids      = todos[i][0]
        orden    = np.argsort(ids)
        todos[i] = todos[i].T[orden].T
        com      = np.in1d(total_id, todos[i][0])

        allcat[:,i*4:i*4+4][com] = todos[i].T

    #Genera el header
    hdr = []
    hdra = 'ID_REF,RA,DEC,X,Y,MAG,MAG_ERR'
    hdr.append(hdra.split(','))
    for i in xrange(len(todos)):
        hdrb = 'ID_1,DX_1,DY_1,NEI_1'
        hdrb = hdrb.replace('1','%d' % nro_epoca[i])
        hdrb = hdrb.split(',')
        hdr.append(hdrb)

    hdr = np.hstack(hdr).tolist()

    allcat = np.hstack([refdata, allcat])
    no_ids = np.isfinite(allcat[:,0])
    allcat = allcat[no_ids]
    nans = np.isnan(allcat)
    allcat[nans] = -9898

    output = Table(allcat, names=hdr)
    print 'Guardando datos...'
    output.write('PM.hdf5', path='data', compression=True)
    #ascii.write(output, 'PM.dat', delimiter=',', fill_values=[('-9898','')])

else:
    print '\nPM.dat encontrado, no se creo archivo!'

sys.exit(1)
#Calcula los PM
yr  = np.genfromtxt('zinfo_img',unpack=True,usecols=(6,))
see = np.genfromtxt('zinfo_img',unpack=True,usecols=(4,))
yep = np.genfromtxt('zinfo_img',unpack=True,usecols=(0,),dtype='string')

yr_ma = np.array(['k' in y for y in yep])
yr    = yr[yr_ma]
see   = see[yr_ma]

yrs = (yr-yr[0])/365.25
yrs = yrs[nro_epoca-1]
see = see[nro_epoca-1]

dxdy_data = np.array(ascii.read('PM.dat', format='csv', fill_values=('',np.nan)))
dxdy_data = np.array(dxdy_data.tolist())

ids = dxdy_data[:,0]
mag = dxdy_data[:,5]
ra  = dxdy_data[:,1]
dec = dxdy_data[:,2]
nei = dxdy_data[:,10::4]
dx  = dxdy_data[:,8::4]
dy  = dxdy_data[:,9::4]

nei_sum  = np.sum(np.isnan(nei), axis=1)
nei_mas  = nei_sum == nei.shape[1]
nei[:,0][nei_mas] = 999

mean_nei = np.nanmean(nei, axis=1)
std_nei  = np.nanstd(nei, axis=1)

#mean_nei[np.isnan(mean_nei)] = 999
#std_nei[np.isnan(std_nei)]   = 999

dx_fin = np.isfinite(dx)*(dx!=888.8)
dy_fin = np.isfinite(dy)*(dy!=888.8)

PM_X = np.zeros(dx_fin.shape[0]) - 999
PM_Y = np.zeros(dy_fin.shape[0]) - 999

count = np.sum(np.isfinite(dx), axis=1)

def PM_calc(i):
    if not dx_fin[i].sum() >= nframes:
        return np.nan, np.nan, np.nan, np.nan
    else:
        ma  = dx_fin[i]
        x   = yrs[ma]
        yx  = dx[i][ma]
        yy  = dy[i][ma]
        ss  = see[ma]

        #Algoritmo antiguo
        '''
        #Ajusta pmx
        popt, pcov = curve_fit(recta, x, yx, sigma=ss)
        pmxx = popt[0]
        pmex = np.sqrt(pcov[0,0])

        #Ajusta pmy
        popt, pcov = curve_fit(recta, x, yy, sigma=ss)
        pmyy = popt[0]
        pmey = np.sqrt(pcov[0,0])
        '''

        #Bayesian Linear Regression (no funciona por ahora y mas costoso)
        '''
        br = BayesianRegression()
        br.fit(x[:, np.newaxis], yx)

        pmxx = br.coef_[0]
        pmex = br.beta_

        br.fit(x[:, np.newaxis], yy)

        pmyy = br.coef_[0]
        pmex = br.beta_
        '''

        #Con matrices
        mu, sig = linear_regression(x, yx, ss)
        pmxx    = mu[0]
        pmex    = sig[0,0]

        mu, sig = linear_regression(x, yy, ss)
        pmyy    = mu[0]
        pmey    = sig[0,0]

        return pmxx, pmyy, pmex, pmey

#PMS_all = np.transpose(Parallel(n_jobs=cpun/2, verbose=8)(delayed(PM_calc)(i) for i in xrange(len(dx))))
PMS_all = np.transpose(barra(PM_calc, xrange(len(dx)),nprocs))

#PMS_all = []
#for i in ProgressBar(xrange(len(dx))):
#    apm = PM_calc(i)
#    PMS_all.append(apm)

#PMS_all = np.transpose(PMS_all)


PMS = np.array(PMS_all[0:2])
PME = np.array(PMS_all[2:])
PMS = PMS * 1000 * 0.339
PME = PME * 1000 * 0.339

PM_X, PM_Y = PMS
PMXE, PMYE = PME
#PM_X, PM_Y = PM_X[np.isfinite(PM_X)], PM_Y[np.isfinite(PM_Y)]
#ids = ids[np.isfinite(PM_X)]

#Plots
pmxa = PM_X[np.isfinite(PM_X)]
pmya = PM_Y[np.isfinite(PM_Y)]

#nbins = int(np.sqrt(len(pmxa)))
nbins = np.arange(-limplot, limplot+nbins, nbins)
#nbins    = np.arange(-30, 30, 0.1)

fig, ax = plt.subplots()
ax.plot(PM_X, PM_Y, '.k', alpha=.75, ms=2)
ax.set(xlim=(-limplot, limplot), ylim=(-limplot, limplot))
ax.text(-25, 25, 'Nro estrellas: %d' % np.isfinite(PM_X).sum())
plt.savefig('VPD.pdf', dpi=200)

H, xedges, yedges = np.histogram2d(pmxa, pmya, bins=nbins)
H  = np.rot90(H)
H  = np.flipud(H)
Hm = np.ma.masked_where(H==0, H)

gs   = gridspec.GridSpec(2, 2, width_ratios=[3,1], height_ratios=[1,3])
figh = plt.figure(figsize=[8,7])
ax2  = plt.subplot(gs[2])
axu  = plt.subplot(gs[0])
axr  = plt.subplot(gs[3])

ax2.pcolormesh(xedges, yedges, Hm, cmap='hot')
axu.hist(pmxa, bins=nbins, histtype='step')
axr.hist(pmya, bins=nbins, histtype='step', orientation='horizontal')

axu.set(xlim=(-limplot, limplot))
axr.set(ylim=(-limplot, limplot))
ax2.set_xlim(-limplot, limplot)
ax2.set_ylim(-limplot, limplot)

plt.savefig('VPDH.pdf', dpi=200)

PM_X[np.isnan(PM_X)] = 999
PM_Y[np.isnan(PM_Y)] = 999
PMXE[np.isnan(PMXE)] = 999
PMYE[np.isnan(PMYE)] = 999

fmt = '%d %.6f %.6f %.6f %.6f %.3f %d %.6f %.6f %.0f %.2f'
hdr = 'ID RA DEC PM_X PM_Y MAG_K NFRAMES PMXE PMYE NEI NEI_STD'

if not os.path.isfile('PM_final.dat'):
    np.savetxt('PM_final.dat', np.transpose([ids, ra, dec, PM_X, PM_Y, mag, count, PMXE, PMYE, mean_nei, std_nei]), fmt=fmt, header=hdr)
else:
    print '\nPM_final.dat encontrado, no se creo archivo!'

#plt.show()

print 'Done!'
