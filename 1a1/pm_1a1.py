import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt

match = False

#Lee los archivos con DX DY
nro_ref    = int(sys.argv[1])
referencia = glob.glob('*k*%03d.dat' % nro_ref)[0]

archivos   = glob.glob('PM_*')
nro_arch   = len(archivos)

#Realiza un match de los IDs de referencia con los DX DY
if match:
    ejecuta  = 'java -jar stilts.jar tmatchn multimode=pairs nin=%d matcher=exact ' % nro_arch
    ejecuta += 'in1=%s ifmt1=ascii values1="ID" ' % referencia

    for i in range(1,nro_arch):
        ejecuta += 'in%d=%s ifmt%d=ascii values%d="ID" ' % (i+1,archivos[i-1],i+1,i+1)

    ejecuta += 'join1=match out=PM.dat ofmt=ascii'

    os.system(ejecuta)

#Calcula los PM
yr = np.genfromtxt('zinfo_img',unpack=True,usecols=(6,),skiprows=6)
yrs = (yr-yr[0])/365.25
yrs = yrs[1:nro_arch]

dxdy_data = np.genfromtxt('PM.dat')

ids = dxdy_data[:,0]
mag = dxdy_data[:,5]
dx  = dxdy_data[:,8::3]
dy  = dxdy_data[:,9::3]

dx_fin = np.isfinite(dx)
dy_fin = np.isfinite(dy)

PM_X = np.zeros(dx_fin.shape[0]) - 999
PM_Y = np.zeros(dy_fin.shape[0]) - 999

for i in range(dx_fin.shape[0]):
    #Calcula PM_X
    if dx_fin[i].sum() > 49:
        ma = dx_fin[i]
        x  = yrs[ma]
        y  = dx[i][ma]

        coeffx  = np.polyfit(x,y,1)
        PM_X[i] = coeffx[0]

    #Calcula PM_Y
    if dy_fin[i].sum() > 49:
        ma = dy_fin[i]
        x  = yrs[ma]
        y  = dy[i][ma]

        coeffy  = np.polyfit(x,y,1)
        PM_Y[i] = coeffy[0]

#Deja como NaN donde no hay PM
PM_X[PM_X == -999] = np.nan
PM_Y[PM_Y == -999] = np.nan

#Plot
fig, ax = plt.subplots()
ax.plot(PM_X,PM_Y,'.')
plt.show()

np.savetxt('PM_final.dat',np.transpose([ids,PM_X,PM_Y,mag]))
