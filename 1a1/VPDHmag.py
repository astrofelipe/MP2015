import matplotlib.pyplot as plt
import numpy as np

#PARAMETROS
lim    = 20
mags   = 11, 14
delta  = 1

#PIPELINE

celdas = mags[1] - mags[0]
sep    = np.arange(mags[0], mags[1]+delta, delta)

ids, ra, dec, pmx, pmy, mag = np.genfromtxt('PM_final.dat', unpack=True)

fig, ax = plt.subplots(nrows=celdas, figsize=[4.2, 4.2*celdas])
fip, ap = plt.subplots(nrows=celdas, figsize=[4.2, 4.2*celdas])
bins = np.arange(-20,21,1)
co   = []

for i in range(len(ax)):
    if i == 0:
        ma = mag < sep[1]
        ax[i].set_title('Ks < %.1f' % sep[1])
        ap[i].set_title('Ks < %.1f' % sep[1])
    elif i== len(ax)-1:
        ax[i].set_title('Ks > %.1f' % sep[-2])
        ap[i].set_title('Ks > %.1f' % sep[-2])
        ma = mag > sep[-2]
    else:
        ax[i].set_title('%.1f > Ks > %.1f' % (sep[i+1], sep[i]))
        ap[i].set_title('%.1f > Ks > %.1f' % (sep[i+1], sep[i]))
        ma = (mag > sep[i])*(mag < sep[i+1])

    x = pmx[ma]
    y = pmy[ma]

    H, xedges, yedges = np.histogram2d(x,y,bins=bins)

    H = np.rot90(H)
    H = np.flipud(H)
    Hmasked = np.ma.masked_where(H==0,H)

    co.append(ax[i].pcolormesh(xedges,yedges,Hmasked, cmap='coolwarm'))
    plt.colorbar(co[i], ax=ax[i])

    ap[i].plot(x, y, '.k', ms=2)

    ax[i].set(xlim=(-lim, lim), ylim=(-lim, lim))
    ap[i].set(xlim=(-lim, lim), ylim=(-lim, lim))

fig.savefig('VPDHmag.pdf', dpi=200, bbox_inches='tight')
fip.savefig('VPDmag.pdf', dpi=200, bbox_inches='tight')
