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

    print('Extrayendo DX y DY')
    idx = np.where(idx)[0][0]
    row = f['data'][idx]
    dx  = row[8::4]
    dy  = row[9::4]

    output = np.vstack((pms, dx, dy)).T
    if hidenan:
        print('No se guardaran filas donde la estrella no esta')
        output = output[np.isfinite(dx)]

    print('Guardando...')
    np.savetxt('id_estrella_de_interes.dat', output, fmt='%03d %.3f %.3f', header='EPOCH DX DY')
