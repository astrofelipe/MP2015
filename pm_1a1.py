from __future__ import division, print_function

import os
import sys
import glob
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import multiprocessing
import subprocess
from matplotlib import gridspec
from scipy.optimize import curve_fit
from astropy.io import ascii
from astropy.utils.console import ProgressBar, color_print
from astropy.table import Table, join, hstack
from astropy.stats import mad_std
from pm_funcs import barra
from linfit import linfit
import pm_funcs

#PARAMETROS
nframes, nbins, limplotpm, nprocs, sig_iter, nsigma, weight = pm_funcs.get_pm1a1()
limplot = limplotpm

color_print('[PM_1a1.py]', 'yellow')
color_print('Parametros:', 'lightgray')
print('\tnframes: %d' % nframes)
print('\tnbins:   %d' % nbins)
print('\tnprocs:  %d' % nprocs)

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

nro_epoca = np.sort([int(f.split('_')[1].split('.')[0]) for f in archivos])
color_print('Epocas:', 'lightgray')
print(nro_epoca)

#Realiza el match entre los PM_*.dat
if not os.path.isfile('PM.h5'):
    print('\nPM.h5 no encontrado, generando archivo...')
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
    #allcat    = np.zeros((maximo, len(todos)*4))
    allcat    = np.memmap('allcat.temp', dtype='float32', mode='w+', shape=(maximo, len(todos)*6))
    allcat[:] = np.nan


    #Rellena la matriz
    print('Ingresando datos a la matriz...')
    for i in ProgressBar(xrange(len(todos))):
        ids      = todos[i][0]
        orden    = np.argsort(ids)
        todos[i] = todos[i].T[orden].T
        com      = np.in1d(total_id, todos[i][0])

        allcat[:,i*6:i*6+6][com] = todos[i].T

    del ids, orden, com

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
    #allcat[nans] = -9898

    #output = Table(allcat, names=hdr)
    print('Guardando PM.h5...')
    h5f = h5py.File('PM.temp', 'w')
    h5f.create_dataset('data', data=allcat)
    h5f.close()
    subprocess.call('mv PM.temp PM.h5', shell=True)
    dxdy_data = allcat
    #output.write('PM.hdf5', path='data', compression=True)

    #print 'Guardando PM.dat...'
    #nans = np.isnan(allcat)
    #allcat[nans] = -9898
    #output = Table(allcat, names=hdr)

    #output.write('PM.dat', fill_values=[('-9898','')], format='ascii.csv')
    #del output, nans, no_ids, hdr, total_id, todos, allcat

else:
    print('\nPM.h5 encontrado, no se creo archivo!')
    print('\nAbriendo PM.h5')
    h5f = h5py.File('PM.h5', 'r')
    dxdy_data = h5f['data']
    #dxdy_data = np.array(Table.read('PM.hdf5', path='data'))
if os.path.isfile('PM_final.dat'):
    print('\nPM_final.dat encontrado! Bye')
    sys.exit(1)


#Calcula los PM

#Lee zinfo
yr  = np.genfromtxt('zinfo_img',unpack=True,usecols=(6,))
see = np.genfromtxt('zinfo_img',unpack=True,usecols=(4,))
yep = np.genfromtxt('zinfo_img',unpack=True,usecols=(0,),dtype='string')

#Filtra solo la banda K
yr_ma = np.array(['k' in y for y in yep])
yep   = yep[yr_ma]
yr    = yr[yr_ma]
see   = see[yr_ma]


#Quiero solo las epocas con las que estoy trabajando
#Busca si el numero de epoca esta en yep
molde = yep[0].split('-')[0]
tengo = np.array([molde+'-%03d.fits' % i for i in nro_epoca])
eff_epoch = np.in1d(yep, tengo) #Epocas que tengo de zinfo
eff_tengo = np.in1d(tengo, yep) #Epocas de zinfo que tengo

yrs = (yr-yr[0])/365.2422 #yr[0] deberia dar igual, siempre que importe solo la pendiente
yrs = yrs[eff_epoch]
see = see[eff_epoch]

#Recupero los NaN de la parte previa para que sea mas facil ignorarlos
#dxdy_data = np.array(dxdy_data.tolist())
#dxdy_data[dxdy_data==-9898] = np.nan

#Identificar columnas
ids = dxdy_data[:,0]
mag = dxdy_data[:,5]
ra  = dxdy_data[:,1]
dec = dxdy_data[:,2]
nei = dxdy_data[:,10::6]
dx  = dxdy_data[:,8::6]
dy  = dxdy_data[:,9::6]
dxe = dxdy_data[:,11::6]
dye = dxdy_data[:,12::6]

#Saca las columnas que no estan en el zinfo (no tendre informacion de yrs)
nei = nei.T[eff_tengo].T
dx  = dx.T[eff_tengo].T
dy  = dy.T[eff_tengo].T
dxe = dxe.T[eff_tengo].T
dye = dye.T[eff_tengo].T

#Obtengo el numero de vecinos usados y pongo 999 los que no cumplen la condicion
nei_sum  = np.sum(np.isnan(nei), axis=1)
nei_mas  = nei_sum == nei.shape[1]
nei[:,0][nei_mas] = 999

mean_nei = np.nanmean(nei, axis=1)
std_nei  = np.nanstd(nei, axis=1)

#Mascara para los dx o dy que tienen valores validos (ie no NaN ni 888)
dx_fin = np.isfinite(dx)*(~np.isclose(dx, 888.8, rtol=1e-4))
dy_fin = np.isfinite(dy)*(~np.isclose(dy, 888.8, rtol=1e-4))

#Aqui se guardan los PM en la direccion X e Y
PM_X = np.zeros(dx_fin.shape[0]) - 999
PM_Y = np.zeros(dy_fin.shape[0]) - 999

count = np.sum(np.isfinite(dx), axis=1)

def PM_calc(i):
    if not dx_fin[i].sum() >= nframes:
        #Si no esta en el minimo de frames, devuelve NaN
        return np.nan, np.nan, np.nan, np.nan, 0
    else:
        ma  = dx_fin[i]*dy_fin[i]
        #print dx.shape, dx_fin.shape, ma.shape
        #print yrs.shape, dx[i].shape, dy[i].shape, see.shape
        x   = yrs[ma]
        yx  = dx[i][ma]
        yxe = dxe[i][ma]
        yy  = dy[i][ma]
        yye = dye[i][ma]
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
        '''
        mu, sig = linear_regression(x, yx, 1.0/ss**2)
        pmxx    = mu[0]
        pmex    = sig[0,0]

        mu, sig = linear_regression(x, yy, 1.0/ss**2)
        pmyy    = mu[0]
        pmey    = sig[0,0]
        '''

        #Con linfit (rapido!)
        if not weight:
            yye = np.ones(len(x))
            yye = np.ones(len(x))

        fitx, cvm = linfit(x, yx, sigmay=yxe)
        pmxx = fitx[0]
        pmex = np.sqrt(cvm[0,0])

        fity, cvm = linfit(x, yy, sigmay=yye)
        pmyy = fity[0]
        pmey = np.sqrt(cvm[0,0])

        ##Sigma Clip
        clip = np.ones(len(x)).astype(bool)
        if sig_iter>0:
            for si in range(sig_iter):
                modelx = recta(x, *fitx)
                resx   = yx - modelx
                stdx   = mad_std(resx[np.isfinite(resx)])

                modely = recta(x, *fity)
                resy   = yy - modely
                stdy   = mad_std(resy[np.isfinite(resy)])

                res = np.sqrt(resx**2 + resy**2)
                std = np.sqrt(stdx**2 + stdy**2)
                clip = clip * (res <= nsigma*std)

                if clip.sum() < nframes:
                    continue

                #Vuelve a calcular
                fitx, cvm = linfit(x[clip], yx[clip], sigmay=yxe[clip])
                pmxx = fitx[0]
                pmex = np.sqrt(cvm[0,0])

                fity, cvm = linfit(x[clip], yy[clip], sigmay=yye[clip])
                pmyy = fity[0]
                pmey = np.sqrt(cvm[0,0])

        return pmxx, pmyy, pmex, pmey, clip.sum()

#PMS_all = np.transpose(Parallel(n_jobs=cpun/2, verbose=8)(delayed(PM_calc)(i) for i in xrange(len(dx))))
color_print('Calculando PMs...', 'orange')
PMS_all = np.transpose(barra(PM_calc, xrange(len(dx)),nprocs))

#Separo PM y errores, ademas escala y convierte a mas/yr
PMS = np.array(PMS_all[0:2])
PME = np.array(PMS_all[2:4])
PMS = PMS * 1000 * 0.339
PME = PME * 1000 * 0.339

clip_epoch = np.array(PMS_all[-1])

PM_X, PM_Y = PMS
PMXE, PMYE = PME

#Plots
pmxa = PM_X[np.isfinite(PM_X)*np.isfinite(PM_Y)]
pmya = PM_Y[np.isfinite(PM_X)*np.isfinite(PM_Y)]

nbins = np.arange(-limplot, limplot+nbins, nbins)

fig, ax = plt.subplots()
ax.plot(PM_X, PM_Y, '.k', alpha=.75, ms=2, rasterized=True)
ax.set(xlim=(-limplot, limplot), ylim=(-limplot, limplot))
ax.text(-15, 15, 'Nro estrellas: %d' % np.isfinite(PM_X).sum())
print('Guardando VPD.png')
plt.savefig('VPD.png', dpi=200)

print(pmxa.shape, pmya.shape)
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

print('Guardando VPDH.png')
plt.savefig('VPDH.png', dpi=200)

PM_X[np.isnan(PM_X)] = 999
PM_Y[np.isnan(PM_Y)] = 999
PMXE[np.isnan(PMXE)] = 999
PMYE[np.isnan(PMYE)] = 999

fmt = '%d %.6f %.6f %.6f %.6f %.3f %d %d %.6f %.6f %.0f %.2f'
hdr = 'ID RA DEC PM_X PM_Y MAG_K NFRAMES CFRAMES PMXE PMYE NEI NEI_STD'

print('Guardando PM_final.dat')
np.savetxt('PM_final.dat', np.transpose([ids, ra, dec, PM_X, PM_Y, mag, count, clip_epoch, PMXE, PMYE, mean_nei, std_nei]), fmt=fmt, header=hdr)

print('Done!')
