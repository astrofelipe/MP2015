import numpy as np
import sys
import os
import pm_funcs

#PARAMETROS
match_tol, modo_ms = pm_funcs.get_master_stilts()

print 'Usando...'
print 'match_tol: %.3f' % match_tol
print 'modo:      %s' % modo_ms

#FIN PARAMETROS
stilts_folder = os.path.dirname(os.path.realpath(__file__)) #STILTS debe estar con el archivo .py

#Lee el archivo con los Inputs
inputs = sys.argv[1]
archivos = np.genfromtxt(inputs, dtype='string')

#Selecciona modo
matcher = 'sky params=%.3f' % match_tol
if modo_ms == 'ID':
    matcher = 'exact'

if modo_ms == 'RA DEC':
    modo_ms = '"RA DEC"'

supermatch  = 'java -jar %s/stilts.jar tmatchn multimode=group nin=%d matcher=%s ' % (stilts_folder, len(archivos), matcher)
for i in xrange(1, len(archivos)+1):
    supermatch += 'in%d=%s ifmt%d=ascii values%d=%s join%d=always ' % (i, archivos[i-1], i, i, modo_ms, i)
supermatch += 'out=master_stilts.dat ofmt=csv'

print supermatch

os.system(supermatch)
