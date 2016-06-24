from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
from astroML.density_estimation import XDGMM
from astroML.plotting.tools import draw_ellipse
#from extreme_deconvolution import extreme_deconvolution

#Parametos (pasar a zparams)
min_mag = 5
max_mag = 25    #Magnitudes para cortar
max_err = 1

compo = 3    #Numero de componentes
itera = 50   #Numero de iteraciones

#Lee
ids, pmx, pmy, magK, pmex, pmey = np.genfromtxt('PM_final.dat', unpack=True, usecols=(0,3,4,5,8,9))

#Filtra
mag_mask = (magK < max_mag) & (magK > min_mag)
err_mask = (pmex**2 + pmey**2)**0.5 < max_err

dataset = np.vstack([pmx[mag_mask*err_mask], pmy[mag_mask*err_mask]]).T
dataerr = np.zeros(dataset.shape + dataset.shape[-1:])
diag    = np.arange(dataset.shape[-1])
dataerr[:, diag, diag] = np.vstack([pmex[mag_mask*err_mask]**2, pmey[mag_mask*err_mask]**2]).T

clf = XDGMM(compo, itera, verbose=True)
clf.fit(dataset, dataerr)
samples = clf.sample(np.sum(mag_mask*err_mask))


print('Centros:\n', clf.mu)
print('Covarianza:\n', clf.V)
print('Alpha:\n', clf.alpha)

#aaa = extreme_deconvolution(dataset, dataerr, np.ones(3)/2, clf.mu, clf.V)
#print(aaa,aaa.xamp_tmp)#, aaa.xmean, aaa.xcovar)

bins = np.arange(-30, 30, 0.5)

fig, ax = plt.subplots(figsize=[4*3,3*2], nrows=2, ncols=3)
ax[0,0].plot(pmx, pmy, '.k', ms=.5, rasterized=True)
ax[0,1].plot(pmx[mag_mask*err_mask], pmy[mag_mask*err_mask], '.k', ms=.5, rasterized=True)
ax[0,2].plot(samples[:,0], samples[:,1], '.k', ms=.5, rasterized=True)

for i in range(clf.n_components):
    draw_ellipse(clf.mu[i], clf.V[i], scales=[2], ax=ax[1,2],
                 ec='k', fc='r', alpha=0.3)

for a in ax[0]:
    a.plot(clf.mu[:,0], clf.mu[:,1], 'xr')

ax[1,0].hist2d(pmx, pmy, cmap='viridis', bins=bins)
ax[1,1].hist2d(pmx[mag_mask*err_mask], pmy[mag_mask*err_mask], cmap='viridis', bins=bins)

for a in np.ravel(ax):
    a.set_xlim(-29,29)
    a.set_ylim(-29,29)
    a.set_aspect('equal')

#fig.tight_layout()
fig.savefig('XD.png', dpi=200, bbox_inches='tight')
