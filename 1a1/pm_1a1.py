import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from sklearn import linear_model
import multiprocessing
from joblib import Parallel, delayed
from matplotlib import gridspec
from scipy.optimize import curve_fit

#PARAMETROS
nframes = 3    #Numero minimo de epocas en que debe estar la estrella
nbins   = 1   #Delta del bin

match = True
stilts_folder = os.path.dirname(os.path.realpath(__file__))
cpun = multiprocessing.cpu_count()

def recta(x,a,b):
    return a*x + b

#Lee los archivos con DX DY
referencia = sys.argv[1]

archivos   = np.sort(glob.glob('./PMs/PM_*'))
nro_arch   = len(archivos)
#threshold  = int(nro_arch * thr_per / 100)

nro_epoca = np.sort([int(f.split('_')[1].split('.')[0]) for f in archivos])
print 'Epocas: ', nro_epoca

#Realiza el match entre los PM_.dat
if not os.path.isfile('PM.dat'):
    ejecuta  = 'java -jar %s/stilts.jar tmatchn multimode=pairs nin=%d matcher=exact ' % (stilts_folder, nro_arch+1)
    ejecuta += 'in1=%s ifmt1=ascii values1=\"ID\" join1=always ' % referencia
    for i in range(1,nro_arch+1):
        ejecuta += 'in%d=%s ifmt%d=ascii values%d=\"ID\" ' % (i+1, archivos[i-1], i+1, i+1)
    ejecuta += 'out=PM.dat ofmt=ascii'

    os.system(ejecuta)
else:
    print '\nPM.dat encontrado, no se creo archivo!'

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

dxdy_data = np.genfromtxt('PM.dat')

ids = dxdy_data[:,0]
mag = dxdy_data[:,5]
ra  = dxdy_data[:,1]
dec = dxdy_data[:,2]
dx  = dxdy_data[:,8::3]
dy  = dxdy_data[:,9::3]

dx_fin = np.isfinite(dx)*(dx!=888.8)
dy_fin = np.isfinite(dy)*(dy!=888.8)

PM_X = np.zeros(dx_fin.shape[0]) - 999
PM_Y = np.zeros(dy_fin.shape[0]) - 999

count = np.sum(dx_fin, axis=1)

def PM_calc(i):
    if not dx_fin[i].sum() >= nframes:
        return np.nan, np.nan
    else:
        ma  = dx_fin[i]
        x   = yrs[ma]
        yx  = dx[i][ma]
        yy  = dy[i][ma]
        ss  = see[ma]

        #model = linear_model.RANSACRegressor(linear_model.LinearRegression())
        #model = linear_model.LinearRegression()
        #gp = GaussianProcess(theta0=1e-4, nugget=1e-10)
        #gp.fit(x[:, np.newaxis], yx)
        #x_pred = np.linspace(x.min(), x.max(),1000)[:, np.newaxis]
        #y_pred, mse = gp.predict(x[:, np.newaxis], eval_MSE=True)
        #gsig = np.sqrt(mse)
        #conf = 2.236 * gsig
        #clip = np.abs(y_pred - yx) < conf

        #Ajusta pmx
        #model.fit(x[:, np.newaxis], yx)
        popt, pcov = curve_fit(recta, x, yx, sigma=ss)
        #pmxx = model.estimator_.coef_[0,0]
        pmxx = popt[0]

        #Ajusta pmy
        #model.fit(x[:, np.newaxis], yy)
        popt, pcov = curve_fit(recta, x, yy, sigma=ss)
        #pmyy = model.estimator_.coef_[0,0]
        pmyy = popt[0]

        return pmxx, pmyy

PMS = np.transpose(Parallel(n_jobs=cpun/2, verbose=8)(delayed(PM_calc)(i) for i in range(len(dx))))
PMS = PMS * 1000 * 0.339
PM_X, PM_Y = PMS
#PM_X, PM_Y = PM_X[np.isfinite(PM_X)], PM_Y[np.isfinite(PM_Y)]
#ids = ids[np.isfinite(PM_X)]

#Plots
pmxa = PM_X[np.isfinite(PM_X)]
pmya = PM_Y[np.isfinite(PM_Y)]

#nbins = int(np.sqrt(len(pmxa)))
nbins = np.arange(-30, 30+nbins, nbins)
#nbins    = np.arange(-30, 30, 0.1)

fig, ax = plt.subplots()
ax.plot(PM_X, PM_Y, '.k', alpha=.75, ms=2)
ax.set(xlim=(-30, 30), ylim=(-30, 30))
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

axu.set(xlim=(-30, 30))
axr.set(ylim=(-30, 30))
ax2.set_xlim(-30, 30)
ax2.set_ylim(-30, 30)

plt.savefig('VPDH.pdf', dpi=200)

PM_X[np.isnan(PM_X)] = 999.9
PM_Y[np.isnan(PM_Y)] = 999.9

fmt = '%d %.6f %.6f %.6f %.6f %.3f %d'
hdr = 'ID RA DEC PM_X PM_Y MAG_K NFRAMES'

if not os.path.isfile('PM_final.dat'):
    np.savetxt('PM_final.dat', np.transpose([ids, ra, dec, PM_X, PM_Y, mag, count]), fmt=fmt, header=hdr)
else:
    print '\nPM_final.dat encontrado, no se creo archivo!'

plt.show()
