import numpy as np
import sys
import os

match_tol     = .3
modo          = 'RA DEC'                                     #Columnas para hacer el match (ej 'RA DEC' o 'ID')
stilts_folder = os.path.dirname(os.path.realpath(__file__)) #STILTS debe estar con el archivo .py

#Lee el archivo con los Inputs
inputs = sys.argv[1]
archivos = np.genfromtxt(inputs, dtype='string')

supermatch  = 'java -jar %s/stilts.jar tmatchn multimode=group nin=%d matcher=sky params=%.3f ' % (stilts_folder, len(archivos), match_tol)
for i in range(1, len(archivos)+1):
    supermatch += 'in%d=%s ifmt%d=ascii values%d=%s join%d=always ' % (i, archivos[i-1], i, modo, i, i)
supermatch += 'out=master_stilts.dat ofmt=csv'

print supermatch

os.system(supermatch)
