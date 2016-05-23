import sys
import h5py

theid = sys.argv[1]

with h5py.File('PM.h5') as f:
    ids = f['data'][:,0]
    print ids
