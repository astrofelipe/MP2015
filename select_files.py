# Genera listas de archivos

# Ejecutar como
# python select_files.py <keywords> <output> [<rango de archivos> <intervalo de archivos>]

# Los parametros entre [] son opcionales, pero deben ir ambos si es que se usan

# El rango de archivos debe tener el formato "archivo_inicial-archivo_final" ej. "4-10" (sin comillas)

# El intervalo de archivos define cada cuantas epocas va a guardar, ej "4-10 2"
# guarda los archivos encontrados numero 4, 6, 8 y 10

import glob
import sys
import numpy as np

keywords = sys.argv[1]
output   = sys.argv[2]

if len(sys.argv) > 3:
    try:
        rango = sys.argv[3]
        salto = sys.argv[4]
    except IndexError:
        print 'Numero de parametros invalidos!'
        sys.exit()

    try:
        salto = int(salto)
        ep1, ep2 = rango.split('-')
        ep1, ep2 = int(ep1), int(ep2)
    except ValueError:
        print 'Parametros invalidos!'
        sys.exit()

archivos = glob.glob(keywords)
print 'Filtrando %d archivos encontrados...' % len(archivos)
archivos = archivos[ep1-1:ep2+1:salto]

print 'Guardando lista con %d archivos...' % len(archivos)
np.savetxt(output,archivos,fmt='%s')
