import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

#Parametros
mag1, mag2 = 11, 14
divisiones = 5#mag2-mag1
limit	   = .05

#Abre datos
pm_file = sys.argv[1]
ids, x, y, col, mag, pmx, pmy = np.genfromtxt(pm_file, unpack=True, usecols=range(7))
ma = np.isfinite(pmx)

ids,x,y,col,mag,pmx,pmy = np.transpose([ids,x,y,col,mag,pmx,pmy])[ma].T

#Crea las divisiones
magdiv = np.linspace(mag1,mag2,divisiones+1)
mask   = [(mag > magdiv[i])*(mag < magdiv[i+1]) for i in range(divisiones)]

#Bines
nbins = 25
bins  = np.linspace(-limit,limit,nbins)

#Plot
fig, ax = plt.subplots(nrows=divisiones, figsize=[2,2*divisiones])
fig2, ax2 = plt.subplots(nrows=divisiones, figsize=[2,2*divisiones])
fig3, ax3 = plt.subplots(nrows=divisiones, figsize=[2,2*divisiones])

fig.subplots_adjust(hspace=0)
fig2.subplots_adjust(hspace=0)
fig3.subplots_adjust(hspace=0)

for i, a in enumerate(ax):
    pmxx = pmx[mask[i]]
    pmyy = pmy[mask[i]]

    a.plot(pmxx,pmyy,'.k',ms=2)
    a.set_xlim(-limit,limit)
    a.set_ylim(-limit,limit)

    H, xedges, yedges = np.histogram2d(pmxx, pmyy, bins=bins)
    H = np.rot90(H)
    H = np.flipud(H)
    Hmasked = np.ma.masked_where(H==0,H)

    ax2[i].pcolormesh(xedges,yedges,Hmasked,cmap='Blues')
    ax2[i].set_xlim(-limit,limit)
    ax2[i].set_ylim(-limit,limit)

    xx, yy = np.mgrid[-limit:limit:100j, -limit:limit:100j]
    pos    = np.vstack([xx.ravel(), yy.ravel()])
    val    = np.vstack([pmxx, pmyy])
    kernel = gaussian_kde(val)
    f      = np.reshape(kernel(pos).T, xx.shape)

    #cfset = ax3[i].contourf(xx, yy, f, 13, cmap='Blues')
    #cset  = ax3[i].contour(xx, yy, f, 13, colors='k', lw=.5)
    #ax3[i].clabel(cset, inline=1)

#fig.tight_layout()
plt.show()
