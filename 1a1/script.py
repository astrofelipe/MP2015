import sys
import os
import subprocess
import numpy as np
from astropy.utils.console import color_print

#PARAMETROS
radio = 3       #Radio en arcsec para iterar
itera = 3       #Numero de iteraciones

print 'Iniciando script con...'
print 'radio: %d' % radio
print 'itera: %d' % itera

#FUNCIONES Y OTROS PARAMETROS
inputs        = sys.argv[1]
ref_cat       = sys.argv[2]
stilts_folder = os.path.dirname(os.path.realpath(__file__))

def makedir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


#PIPELINE
color_print('\nIniciando script, recuerda eliminar la carpeta PMs, el archivo PM.dat, etc. antes de ejecutar','yellow')

output = subprocess.check_output('grep "PDF de Output" %s/tlineal_1a1.py' % stilts_folder, shell=True)
output = output.split("'")[1]

refstars = subprocess.check_output('grep "Catalogo con las estrellas de referencia" %s/tlineal_1a1.py' % stilts_folder, shell=True)
refstars = refstars.split("'")[1]

nframes = subprocess.check_output('grep "Numero minimo de epocas en que debe estar la estrella" %s/pm_1a1.py' % stilts_folder, shell=True)
nframes = int(nframes.split(' ')[2])

print '\nOutput en tlineal_1a1.py: %s' % output
print 'Archivo de refstars: %s' % refstars
print 'nframes en pm_1a1.py: %d' % nframes

if not os.path.isfile('refstars0.gc'):
    print '\nArchivo refstars0.gc no encontrado!'
    sys.exit()

subprocess.call('cp refstars0.gc refstars.gc', shell=True)

for i in range(itera):
    color_print('\nComenzando iteracion: %d' % (i+1), 'lightcyan')

    #Crea carpeta para guardar los outputs
    makedir('iter_%d' % (i+1))

    color_print('\tEjecutando tlineal_1a1.py', 'cyan')
    subprocess.call('python %s/tlineal_1a1.py %s' % (stilts_folder, inputs), shell=True)

    color_print('\tEjecutando pm_1a1.py', 'cyan')
    subprocess.call('python %s/pm_1a1.py %s' % (stilts_folder, ref_cat), shell=True)

    color_print('\tEjecutando VPDHmag.py', 'cyan')
    subprocess.call('python %s/VPDHmag.py' % stilts_folder, shell=True)

    color_print('\tMoviendo archivos', 'cyan')
    subprocess.call('mv %s*.pdf iter_%d' % (output, (i+1)), shell=True)
    subprocess.call('mv %s iter_%d' % (refstars, (i+1)), shell=True)
    makedir('iter_%d/PMs' % (i+1))
    subprocess.call('mv PM*.dat iter_%d' % (i+1), shell=True)
    subprocess.call('mv PMs/* iter_%d/PMs' % (i+1), shell=True)
    subprocess.call('mv VPD*.pdf iter_%d' % (i+1), shell=True)

    color_print('\tGenerando nuevo archivo de refstars', 'cyan')
    ids, pmx, pmy, nf = np.genfromtxt('iter_%d/PM_final.dat' % (i+1), unpack=True, usecols=(0,3,4,6))
    pmr = np.sqrt(pmx**2 + pmy**2)

    mask = (pmr <= radio) * (nf >= nframes)

    data = np.genfromtxt('iter_%d/PM_final.dat' % (i+1))
    data = data[mask]

    fmt = '%d %.6f %.6f %.6f %.6f %.3f %d %.6f %.6f'
    hdr = 'ID RA DEC PM_X PM_Y MAG_K NFRAMES PMXE PMYE'

    np.savetxt(refstars, data, fmt=fmt, header=hdr)

subprocess.call('rm -r PMs', shell=True)
