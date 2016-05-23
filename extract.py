from __future__ import division, print_function
import sys
import h5py
import numpy as np

theid = int(sys.argv[1])

with h5py.File('PM.h5') as f:
    ids = f['data'][:,0]
    idx = ids == theid

    if np.sum(idx) == 0:
        print('ID no encontrado!')
        sys.exit(1)

    idx = np.where(idx)[0][0]
    row = f['data'][idx]
    print(row[:15])
