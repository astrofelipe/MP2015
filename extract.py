from __future__ import division, print_function
import sys
import h5py
import glob
import numpy as np

theid = int(sys.argv[1])
hidenan = True

print('Abriendo archivo...')
with h5py.File('PM.h5') as f:
    print('Buscando ID')
    ids = f['data'][:,0]
    idx = ids == theid

    if np.sum(idx) == 0:
        print('ID no encontrado!')
        sys.exit(1)

    print('Buscando archivos de epocas')
    pms = glob.glob('PMs/*.dat')
    pms = np.sort([int(s.split('_')[-1].split('.')[0]) for s in pms])

    print('Extrayendo seeing')
    yep = np.genfromtxt('zinfo_img',unpack=True,usecols=(0,),dtype='string')
    hak = np.array(['k' in y for y in yep])
    yep = yep[hak]

    molde = yep[0].split('-')[0]
    tengo = np.array([molde+'-%03d.fits' % i for i in pms])
    effep = np.in1d(yep, tengo)

    zinfodat = np.genfromtxt('zinfo_img',usecols=(4,6))[hak][effep]
    see, jd  = zinfodat.T

    print('Extrayendo DX y DY')
    idx = np.where(idx)[0][0]
    row = f['data'][idx]
    dx  = row[8::4]
    dy  = row[9::4]

    output = np.vstack((pms, jd, dx, dy, see)).T
    if hidenan:
        print('No se guardaran filas donde la estrella no esta')
        output = output[np.isfinite(dx)]

    print('Guardando...')
    np.savetxt('%d.dat' % theid, output, fmt='%03d %.3f %.3f %.8f', header='EPOCH DX DY SEEING')
