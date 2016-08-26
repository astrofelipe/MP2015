from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
from matplotlib import rc, cm, mlab
from matplotlib.colors import LinearSegmentedColormap
from sklearn import mixture
from gmm import GMM
from astroML.density_estimation import XDGMM
from astroML.plotting.tools import draw_ellipse
#from extreme_deconvolution import extreme_deconvolution

#Parametos (pasar a zparams)
min_mag = 5
max_mag = 25    #Magnitudes para cortar
max_err = 1

compo = 3    #Numero de componentes
itera = 1   #Numero de iteraciones

#Funciones
def gaussian(x, amp, mu, sig, y0):
    return amp * np.exp(-(x-mu)**2 / (2*sig**2)) + y0

def three_gaussian(x, amp1, mu1, sig1,
                      amp2, mu2, sig2,
                      amp3, mu3, sig3, y0):
    return (gaussian(x, amp1, mu1, sig1, 0) +
    gaussian(x, amp2, mu2, sig2, 0) +
    gaussian(x, amp3, mu3, sig3, 0)) + y0


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

clfu = mixture.VBGMM(compo, covariance_type='full', tol=1e-5, n_iter=1000)
clfu.fit((pmx[mag_mask*err_mask])[:,np.newaxis])
meu  = np.hstack(clfu.means_)
stu  = np.hstack(clfu.precs_)[0]
weu  = np.hstack(clfu.weights_)

clfd = XDGMM(compo, itera, verbose=True)
clfd.fit(pmy[mag_mask*err_mask][:,np.newaxis], (pmey[mag_mask*err_mask]**2)[:,np.newaxis,np.newaxis])
samd = clfd.sample(np.sum(mag_mask*err_mask))

print('Centros:\n', clf.mu)
print('Covarianza:\n', clf.V)
print('Alpha:\n', clf.alpha)

bins = np.arange(-30, 30, 0.5)
xx   = np.linspace(-30, 30, 1e3)
cmap  = cm.get_cmap('viridis')
color = cmap(np.linspace(0.33, 1, cmap.N // 2))
cmap  = LinearSegmentedColormap.from_list('Eeeee', color)

fig = plt.figure(figsize=[8,8])
gs  = gridspec.GridSpec(4,4)

axh = plt.subplot(gs[1:,:-1])
axu = plt.subplot(gs[0,:-1])
axd = plt.subplot(gs[1:,-1])

#axh.set_aspect('equal')

#fig, ax = plt.subplots(figsize=[4*2,4*2], nrows=2, ncols=3)
axh.plot(pmx[mag_mask*err_mask], pmy[mag_mask*err_mask], '.k', ms=.5, rasterized=True)
#ax[0,1].plot(pmx[mag_mask*err_mask], pmy[mag_mask*err_mask], '.k', ms=.5, rasterized=True)
#ax[0,2].plot(samples[:,0], samples[:,1], '.k', ms=.5, rasterized=True)

yh, _, _ = axu.hist(pmx[mag_mask*err_mask], histtype='stepfilled', bins=bins, lw=2, alpha=.9, facecolor=cmap(0.1), edgecolor=cmap(0.1))
axu.scatter(bins[:-1]+0.5, yh)
popt, pcov = curve_fit(three_gaussian, yh, bins[:-1]+0.5, p0=[300, -7, 2, 800, 0, 3, 800, 0, 3, 0])

axu.plot(xx, three_gaussian(xx, *popt), color=cmap(0.9))

#for i in range(3):
#    yy = mlab.normpdf(xx, meu[i], np.sqrt(stu[i]))
#    print(yy.max(), pmx[mag_mask*err_mask].max())
#    axu.plot(xx, yy, color=cmap(0.9))
#axu.hist(samu, bins=bins, histtype='step', lw=1, edgecolor=cmap(0.9))
#for y in clfu.mu:
#    axu.axvline(y, color=cmap(0.9))

axd.hist(pmy[mag_mask*err_mask], histtype='stepfilled', bins=bins, lw=2, alpha=.9, facecolor=cmap(0.1), edgecolor=cmap(0.1), orientation='horizontal')
axd.hist(samd, bins=bins, histtype='step', lw=1, edgecolor=cmap(0.9), orientation='horizontal')
for y in clfd.mu:
    axd.axhline(y, color=cmap(0.9))

#for i in range(clf.n_components):
#    draw_ellipse(clf.mu[i], clf.V[i], scales=[2], ax=ax[1,2],
#                 ec='k', fc='r', alpha=0.3)

#for a in ax[0]:
#    a.plot(clf.mu[:,0], clf.mu[:,1], 'xr')

#ax[1,0].hist2d(pmx, pmy, cmap='viridis', bins=bins)
#ax[1,1].hist2d(pmx[mag_mask*err_mask], pmy[mag_mask*err_mask], cmap='viridis', bins=bins)

axh.set_xlim(-29, 29)
axh.set_ylim(-29, 29)
#for a in np.ravel(ax):
#    a.set_xlim(-29,29)
#    a.set_ylim(-29,29)
#    a.set_aspect('equal')

#fig.tight_layout()
fig.savefig('XD.png', dpi=200, bbox_inches='tight')
