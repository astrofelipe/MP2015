import numpy as np
import sys
import os

#PARAMETROS A EDITAR
match_tol     = .3           #Tolerancia para el modo a usar (ID no usa este parametro)
#modo          = '"RA DEC"'    #Columnas para hacer el match (ej 'RA DEC' o 'ID')
modo          = 'ID'         #Columnas para hacer el match (ej 'RA DEC' o 'ID')

#FIN PARAMETROS
stilts_folder = os.path.dirname(os.path.realpath(__file__)) #STILTS debe estar con el archivo .py

#Lee el archivo con los Inputs
inputs = sys.argv[1]
archivos = np.genfromtxt(inputs, dtype='string')

#Selecciona modo
matcher = 'sky params=%.3f' % match_tol
if modo == 'ID':
    matcher = 'exact'

supermatch  = 'java -jar %s/stilts.jar tmatchn multimode=group nin=%d matcher=%s ' % (stilts_folder, len(archivos), matcher)
for i in range(1, len(archivos)+1):
    supermatch += 'in%d=%s ifmt%d=ascii values%d=%s join%d=always ' % (i, archivos[i-1], i, i, modo, i)
supermatch += 'out=master_stilts.dat ofmt=csv'

print supermatch

os.system(supermatch)
