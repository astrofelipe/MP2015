import numpy as np
import sys

nro = int(sys.argv[1])

zinfo  = np.genfromtxt('zinfo_img', unpack=True)
seeing = zinfo[4]
orden  = np.argsort(seeing)
files  = np.genfromtxt('zinfo_img', dtype='string', usecols=(0,))

output = files[orden][:nro]
output = np.array([o.replace('.fits','.datgc') for o in output])
np.savetxt('mfinput', output, fmt='%s')
