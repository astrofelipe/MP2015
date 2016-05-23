from __future__ import division, print_function
import sys
import h5py
import glob
import numpy as np

theid = int(sys.argv[1])

with h5py.File('PM.h5') as f:
    ids = f['data'][:,0]
    idx = ids == theid

    if np.sum(idx) == 0:
        print('ID no encontrado!')
        sys.exit(1)

    pms = glob.glob('PMs/*.dat')
    pms = np.sort([s.split('_')[-1].splt('.')[0] for s in pms])
    print(pms)

    idx = np.where(idx)[0][0]
    row = f['data'][idx]
    dx  = row[8::4]
    dy  = row[9::4]
    print(dx,dy)
