import sys
import os
import subprocess
import numpy as np

#PARAMETROS
radio = 3       #Radio en arcsec para iterar
itera = 3       #Numero de iteraciones

#FUNCIONES Y OTROS PARAMETROS
inputs        = sys.argv[1]
ref_cat       = sys.argv[2]
stilts_folder = os.path.dirname(os.path.realpath(__file__))

def makedir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


#PIPELINE
print '\nIniciando script, recuerda eliminar la carpeta PMs, el archivo PM.dat, etc. antes de ejecutar'

output = subprocess.check_output('grep "PDF de Output" %s/tlineal_1a1.py' % stilts_folder, shell=True)
output = output.split("'")[1]

refstars = subprocess.check_output('grep "Catalogo con las estrellas de referencia" %s/tlineal_1a1.py' % stilts_folder, shell=True)
refstars = refstars.split("'")[1]

nframes = subprocess.check_output('grep "Numero minimo de epocas en que debe estar la estrella" %s/pm_1a1.py' % stilts_folder, shell=True)
nframes = int(nframes.split(' ')[2])

print '\nOutput en tlineal_1a1.py: %s' % output
print 'Archivo de refstars: %s' % refstars
print 'nframes en pm_1a1.py: %d' % nframes

for i in range(itera):
    print '\nComenzando iteracion: %d' % (i+1)

    #Crea carpeta para guardar los outputs
    makedir('iter_%d' % (i+1))

    print '\tEjecutando tlineal_1a1.py'
    subprocess.call('python %s/tlineal_1a1.py %s' % (stilts_folder, inputs), shell=True)

    print '\tEjecutando pm_1a1.py'
    subprocess.call('python %s/pm_1a1.py %s' % (stilts_folder, ref_cat), shell=True)

    print '\tMoviendo archivos'
    subprocess.call('mv %s*.pdf iter_%d' % (output, (i+1)), shell=True)
    subprocess.call('mv %s iter_%d' % (refstars, (i+1)), shell=True)
    subprocess.call('mv PM*.dat iter_%d' % (i+1), shell=True)
    subprocess.call('mv PMs iter_%d' % (i+1), shell=True)

    print '\tGenerando nuevo archivo de refstars'
    ids, pmx, pmy, nf = np.genfromtxt('iter_%d/PM_final.dat' % (i+1), unpack=True, usecols=(0,3,4,6))
    pmr = np.sqrt(pmx**2 + pmy**2)

    mask = (pmr <= radio) * (nf >= nframes)

    data = np.genfromtxt('iter_%d/PM_final.dat' % (i+1))
    data = data[mask]

    fmt = '%d %.6f %.6f %.6f %.6f %.3f %d'
    hdr = 'ID RA DEC PM_X PM_Y MAG_K NFRAMES'

    np.savetxt(refstars, data, fmt=fmt, hdr=hdr)
