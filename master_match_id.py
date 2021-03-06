import numpy as np
import sys
import h5py
import subprocess
from astropy.io import ascii
from astropy.table import Table
from astropy.utils.console import ProgressBar
#from joblib import Parallel, delayed
import pm_funcs
from pm_funcs import barra

if len(sys.argv)==1:
    print 'Para ejecutar'
    print 'python master_match_id.py <archivo con inputs>'
    print 'El match se hara usando la columna ID'
    sys.exit()

#Parametros
min_epochs, nprocs = pm_funcs.get_master_match()

#Argumentos
input_list = sys.argv[1]

print 'Iniciando con...'
print 'min_epochs: %d' % min_epochs

archivos   = np.genfromtxt(input_list, dtype='string')

def load_file(archivo):
    data = np.genfromtxt(archivo, unpack=True)
    return data

#todos = []
#for f in archivos:
#    data = np.genfromtxt(f, unpack=True)
#    todos.append(data)

print '\nLeyendo archivos...'
#todos = Parallel(n_jobs=4, verbose=8)(delayed(load_file)(f) for f in archivos)
todos = barra(load_file, archivos, nprocs)

#Encuentra el ID maximo
print 'Encontrando maximo de estrellas'
maximos = np.zeros(len(todos))
for i in xrange(len(todos)):
    maximos[i] = np.max(todos[i][0])
maximo = int(np.max(maximos))
print'\tID Maximo: %d' %maximo

#Genera matriz donde van todos los catalogos
allcat    = np.zeros((maximo, len(todos)*7))
allcat[:] = np.nan

total_id  = np.arange(1,maximo+1)

#Rellena la matriz
print 'Ingresando datos a la matriz...'
for i in ProgressBar(xrange(len(todos))):
    ids = todos[i][0]
    com = np.in1d(total_id, ids)
    orden = np.argsort(ids)

    allcat[:,i*7:i*7+7][com] = todos[i].T[orden]

#Genera el header
hdr = []
for i in xrange(len(todos)):
    hdrb = 'ID_1,RA_1,DEC_1,X_1,Y_1,MAG_1,MAG_ERR_1'
    hdrb = hdrb.replace('1','%d' % (i+1))
    hdrb = hdrb.split(',')
    hdr.append(hdrb)

hdr = np.hstack(hdr).tolist()

#Achica el catalogo
ids = allcat[:, 0::7]

found  = np.sum(np.isfinite(ids), axis=1)
rej    = found >= min_epochs
allcat = allcat[rej]

#print 'Convirtiendo a tabla...'
#nans = np.isnan(allcat)
#allcat[nans] = -9898
print 'Guardando...'
h5f = h5py.File('master_match.temp', 'w')
h5f.create_dataset('data', data=allcat)
h5f.close()
subprocess.call('mv master_match.temp master_match.h5', shell=True)
#output = Table(allcat, names=hdr)
#ascii.write(output, 'master_match.dat', delimiter=',', fill_values=[('-9898','')])
#output.write('master_match.h5', path='data')
