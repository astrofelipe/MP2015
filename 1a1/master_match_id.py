import numpy as np
import sys
from astropy.io import ascii
from astropy.table import Table

if len(sys.argv)==1:
    print 'Para ejecutar'
    print 'python master_match_id.py <archivo con inputs>'
    print 'El match se hara usando la columna ID'
    sys.exit()

#Argumentos
input_list = sys.argv[1]
archivos   = np.genfromtxt(input_list, dtype='string')

todos = []
for f in archivos:
    data = np.genfromtxt(f, unpack=True)
    todos.append(data)

#Encuentra el ID maximo
maximos = np.zeros(len(todos))
for i in range(len(todos)):
    maximos[i] = np.max(todos[i][0])
maximo = np.max(maximos)

#Genera matriz donde van todos los catalogos
allcat    = np.zeros((maximo, len(todos)*7))
allcat[:] = np.nan

total_id  = np.arange(1,maximo+1)

#Rellena la matriz
for i in range(len(todos)):
    ids = todos[i][0]
    com = np.in1d(total_id, ids)
    orden = np.argsort(ids)

    print todos[i].T.shape

    allcat[:,i*7:i*7+7][com] = todos[i].T[orden]

#Genera el header
hdr = []
for i in range(len(todos)):
    hdrb = 'ID_1,RA_1,DEC_1,X_1,Y_1,MAG_1,MAG_ERR_1'
    hdrb = hdrb.replace('1','%d' % (i+1))
    hdrb = hdrb.split(',')
    hdr.append(hdrb)

hdr = np.hstack(hdr).tolist()

nans = np.isnan(allcat)
allcat[nans] = -9898

output = Table(allcat, names=hdr)
print output
#ascii.write(allcat, 'master_match.dat', delimiter=',', names=hdr, fill_values=[('-9898','')])
ascii.write(output, 'master_match.dat', delimiter=',', fill_values=[('-9898','')])
