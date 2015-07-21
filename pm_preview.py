import sys
import numpy as np
import matplotlib.pyplot as plt

#Parametros
mag1, mag2 = 11, 17
divisiones = mag2-mag1

#Abre datos
pm_file = sys.argv[1]
ids, x, y, col, mag, pmx, pmy = np.genfromtxt(pm_file, unpack=True, usecols=range(7))

#Crea las divisiones
magdiv = np.linspace(mag1,mag2,divisiones)
masks  = [(mag > magdiv[i])*(mag < magdiv[i+1]) for i in range(divisiones - 1)]

print len(masks)
#Plot
fig, ax = plt.subplots(nrows=divisiones, figsize=[4,3*divisiones])
print len(ax)
for i, a in enumerate(ax):
    pmxx = pmx[mask[i]]
    pmyy = pmy[mask[i]]

    a.plot(pmxx,pmyy,'.k')
plt.show()
