import sys
import numpy as np
import matplotlib.pyplot as plt

#Parametros
mag1, mag2 = 11, 15
divisiones = mag2-mag1
limit      = 30

#Abre datos
pm_file = sys.argv[1]
ids, x, y, col, mag, pmx, pmy = np.genfromtxt(pm_file, unpack=True, usecols=range(7))

#Crea las divisiones
magdiv = np.linspace(mag1,mag2,divisiones+1)
mask   = [(mag > magdiv[i])*(mag < magdiv[i+1]) for i in range(divisiones)]

#Bines
nbins = 25
bins  = np.linspace(-limit,limit,nbins)

#Plot
fig, ax = plt.subplots(nrows=divisiones, figsize=[2,2*divisiones])
fig2, ax2 = plt.subplots(nrows=divisiones, figsize=[2,2*divisiones])
plt.subplots_adjust(hspace=0)

for i, a in enumerate(ax):
    pmxx = pmx[mask[i]]
    pmyy = pmy[mask[i]]

    a.plot(pmxx,pmyy,'.k',ms=1)
    a.set_xlim(-limit,limit)
    a.set_ylim(-limit,limit)

    H, xedges, yedges = np.histogram2d(pmxx, pmyy, bins=bins)
    H = np.rot90(H)
    H = np.flipud(H)
    Hmasked = np.ma.masked_where(H==0,H)

    ax2[i].pcolormesh(xedges,yedges,Hmasked)
    ax2[i].set_xlim(-limit,limit)
    ax2[i].set_ylim(-limit,limit)

#fig.tight_layout()
plt.show()
