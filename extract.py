import sys
import h5py
import numpy as np

theid = sys.argv[1]

with h5py.File('PM.h5') as f:
    ids = f['data'][:,0]
    print type(ids)
    idx = ids == theid

    if np.sum(idx) == 0:
        print 'ID no encontrado!'
        sys.exit(1)
