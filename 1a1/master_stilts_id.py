import numpy as np
import sys

if len(sys.argv)==1:
    print 'Para ejecutar'
    print 'python master_stilts_id.py <archivo con inputs>'
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

    allcat[:,i*7:i*7+7][com] = todos[i].T

#Genera el header
hdrb = 'ID_1,RA_1,DEC_1,X_1,Y_1,MAG_1,MAG_ERR_1,'
hdr  = hdrb
for i in range(1,len(todos)):
    hdr += hdrb.replace('1', '%d' % (i+1))
hdr = hdr[:-1]

np.savetxt('master_stilts.dat', allcat, fmt='%f', header=hdr, delimiter=',')
