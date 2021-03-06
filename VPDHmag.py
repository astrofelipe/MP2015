import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pm_funcs
from astropy.stats import mad_std

#PARAMETROS

limplotvp, mags, delta, min_ep, min_nei, sigma_err = pm_funcs.get_VPDHmag()
limplot = limplotvp

print 'Usando...'
print 'limplot:    %d' % limplot
print 'mags:   %.3f,%.3f' % (mags[0], mags[1])
print 'delta:  %.3f' % delta
print 'min_ep: %d' % min_ep

#PIPELINE

celdas = int((mags[1] - mags[0]) / delta)
sep    = np.arange(mags[0], mags[1]+delta, delta)

data = np.genfromtxt('PM_final.dat', unpack=True, usecols=range(11))
#Rejection por 999
ma_999 = data[3] != 999
data = data.T[ma_999].T

#Rejection por vecinos
ma_nei = data[10] >= min_nei

#Rejection por errores
pme = np.sqrt(data[8]**2 + data[9]**2)
ma_pme = np.abs(pme - np.nanmedian(pme)) < sigma_err*mad_std(pme - np.nanmedian(pme))

#Rejection por numero de epocas
ma_nep = data[6] >= min_ep

ids, ra, dec, pmx, pmy, mag, nep, pmxe, pmye, nei, neistd = np.transpose(data.T[ma_nei*ma_pme*ma_nep])

#Rejection por vecinos
ma_nei = nei >= min_nei

fig, ax = plt.subplots(nrows=celdas, figsize=[1.25*4, 1.25*3*celdas])
fip, ap = plt.subplots(nrows=celdas, figsize=[1.25*4, 1.25*4*celdas])
bins = np.arange(-limplot,limplot+1,1)
co   = []

for i in xrange(len(ax)):
    if i == 0:
        ma = (mag > mags[0])*(mag < sep[1])*(pmx != 999.9)
        ax[i].set_title('%.1f > Ks > %.1f' % (sep[1], mags[0]))
        ap[i].set_title('%.1f > Ks > %.1f' % (sep[1], mags[0]))
    elif i== len(ax)-1:
        ax[i].set_title('%.1f > Ks > %.1f' % (sep[-2], mags[1]))
        ap[i].set_title('%.1f > Ks > %.1f' % (sep[-2], mags[1]))
        ma = (mag > sep[-2])*(mag < mags[1])*(pmx != 999.9)
    else:
        ax[i].set_title('%.1f > Ks > %.1f' % (sep[i+1], sep[i]))
        ap[i].set_title('%.1f > Ks > %.1f' % (sep[i+1], sep[i]))
        ma = (mag > sep[i])*(mag < sep[i+1])*(pmx != 999.9)

    x = pmx[ma]
    y = pmy[ma]

    H, xedges, yedges = np.histogram2d(x,y,bins=bins)

    H = np.rot90(H)
    H = np.flipud(H)
    Hmasked = np.ma.masked_where(H==0,H)

    co.append(ax[i].pcolormesh(xedges,yedges,Hmasked, cmap='coolwarm'))
    plt.colorbar(co[i], ax=ax[i])

    ap[i].plot(x, y, '.k', ms=2, rasterized=True)

    nstars = len(x)
    ax[i].text(.066,.933, 'Nro stars: %d' % nstars, transform = ax[i].transAxes)
    ap[i].text(.066,.933, 'Nro stars: %d' % nstars, transform = ap[i].transAxes)

    ax[i].set(xlim=(-limplot, limplot), ylim=(-limplot, limplot))
    ap[i].set(xlim=(-limplot, limplot), ylim=(-limplot, limplot))

fig.tight_layout()
fip.tight_layout()
fig.savefig('VPDHmag.png', dpi=200)
fip.savefig('VPDmag.png', dpi=200)
