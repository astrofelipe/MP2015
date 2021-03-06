import sys
import argparse
import os
import subprocess
import numpy as np
import pm_funcs
from astropy.utils.console import color_print

#PARAMETROS
radio, itera, output, refer, nframes, min_ep, max_err = pm_funcs.get_script()

print 'Iniciando script con...'
print 'radio: %d' % radio
print 'itera: %d' % itera
print

#ARGUMENTOS
parser = argparse.ArgumentParser(description='Script PM VVV')
parser.add_argument('<Input List>', help='Lista con inputs (para tlineal)')
parser.add_argument('<Ref Catalog>', help='Catalogo de referencia (usado por PM_1a1)')
parser.add_argument('-C', action='store_true', help='Continua desde la ultima iteracion')

args = parser.parse_args()

#FUNCIONES Y OTROS PARAMETROS
inputs        = vars(args)['<Input List>']#sys.argv[1]
ref_cat       = vars(args)['<Ref Catalog>']#sys.argv[2]
stilts_folder = os.path.dirname(os.path.realpath(__file__))

#continua = False
#if (len(sys.argv) > 3) and (sys.argv[3]=='-C'):
#    continua = True
continua = args.C
print continua

def makedir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


#PIPELINE
#output = subprocess.check_output('grep "PDF de Output" %s/tlineal_1a1.py' % stilts_folder, shell=True)
#output = output.split("'")[1]

#refstars = subprocess.check_output('grep "Catalogo con las estrellas de referencia" %s/tlineal_1a1.py' % stilts_folder, shell=True)
#refstars = refstars.split("'")[1]

#nframes = subprocess.check_output('grep "Numero minimo de epocas en que debe estar la estrella" %s/pm_1a1.py' % stilts_folder, shell=True)
#nframes = int(nframes.split(' ')[2])

#min_ep = subprocess.check_output('grep "min_ep =" %s/VPDHmag.py' % stilts_folder, shell=True)
#min_ep = int(min_ep.split(' ')[-1])
refstars = refer
print '\nOutput en tlineal_1a1.py: %s' % output
print 'Archivo de refstars: %s' % refer
print 'nframes en pm_1a1.py: %d' % nframes
print 'min_ep en VPDHmag.py: %d' % min_ep
if continua:
    print 'Continuando desde la ultima iteracion'

nro_files = np.genfromtxt(inputs, unpack=True, dtype='string')
if not nro_files.shape:
    nro_files = np.atleast_1d(nro_files)
nro_files = nro_files.size

if nro_files < nframes:
    print '\nNumero de archivos en el input es menor que nframes!'
    sys.exit()

if nro_files < min_ep:
    print '\nNumero de archivos en el input es menor que min_ep!'
    sys.exit()

if not os.path.isfile('refstars0.gc'):
    print '\nArchivo refstars0.gc no encontrado!'
    sys.exit()

if not continua:
    #subprocess.call('cp refstars0.gc refstars.gc', shell=True)
    ids = np.genfromtxt('refstars0.gc', unpack=True, usecols=(0,))

    if os.path.isfile('zelimchisha'):
        print '\nEliminando estrellas de zelimchisha'
        rej_ids = np.genfromtxt('zelimchisha', unpack=True, usecols=(0,))
        id_mask = ~np.in1d(ids, rej_ids)
        print '\tEncontradas %d estrellas' % (~id_mask).sum()
    else:
        id_mask = np.ones(len(ids)).astype(bool)

    data = np.genfromtxt('refstars0.gc')
    data = data[id_mask]

    print 'Guardando refstars para primera iteracion'
    np.savetxt(refstars, data, fmt='%f')

    last_idx = 0

    for i in range(1,itera+1):
        if os.path.exists('iter_%d' % i):
            subprocess.call('rm -r iter_%d' % i, shell=True)

if continua:
    import glob
    iter_dirs = np.sort(glob.glob('iter_*'))

    if len(iter_dirs) == 0:
        print 'No hay carpetas iter!'
        sys.exit(1)
    last_idx  = int(iter_dirs[-1].split('_')[-1])

    if not os.path.isfile('iter_%d/PM_final.dat' % last_idx):
        print '\nNo se encontro PM_final.dat en iter_%d!' % last_idx
        sys.exit(1)

    print '\nGenerando nuevo archivo de refstars a partir de la ultima iteracion'
    ids, pmx, pmy, nf, pmxe, pmye = np.genfromtxt('iter_%d/PM_final.dat' % (last_idx), unpack=True, usecols=(0,3,4,6,8,9))
    pmr = np.sqrt(pmx**2 + pmy**2)
    pme = np.sqrt(np.square(pmxe) + np.square(pmye))

    if os.path.isfile('zelimchisha'):
        rej_ids = np.genfromtxt('zelimchisha', unpack=True, usecols=(0,))
        id_mask = ~np.in1d(ids, rej_ids)
    else:
        id_mask = np.ones(len(ids)).astype(bool)

    mask = (pmr <= radio) * (nf >= nframes) * (pme <= max_err) * (id_mask)

    data = np.genfromtxt('iter_%d/PM_final.dat' % (last_idx))
    data = data[mask]

    fmt = '%d %.6f %.6f %.6f %.6f %.3f %d %d %.6f %.6f %.0f %.2f'
    hdr = 'ID RA DEC PM_X PM_Y MAG_K NFRAMES CFRAMES PMXE PMYE NEI NEI_STD'

    np.savetxt(refstars, data, fmt=fmt, header=hdr)

for i in range(itera):
    color_print('\nComenzando iteracion: %d' % (i+1), 'lightcyan')

    #Crea carpeta para guardar los outputs
    makedir('iter_%d' % (last_idx+i+1))

    color_print('\tEjecutando tlineal_1a1.py', 'cyan')
    subprocess.call('python -u %s/tlineal_1a1.py %s' % (stilts_folder, inputs), shell=True)

    color_print('\tEjecutando pm_1a1.py', 'cyan')
    subprocess.call('python %s/pm_1a1.py %s' % (stilts_folder, ref_cat), shell=True)

    color_print('\tEjecutando VPDHmag.py', 'cyan')
    subprocess.call('python %s/VPDHmag.py' % stilts_folder, shell=True)

    color_print('\tMoviendo archivos', 'cyan')
    subprocess.call('mv %s*.png iter_%d' % (output, (last_idx+i+1)), shell=True)
    subprocess.call('mv %s iter_%d' % (refstars, (last_idx+i+1)), shell=True)
    makedir('iter_%d/PMs' % (last_idx+i+1))
    subprocess.call('mv PM*.h5 iter_%d' % (last_idx+i+1), shell=True)
    subprocess.call('mv PM_final.dat iter_%d' % (last_idx+i+1), shell=True)
    subprocess.call('mv PMs/* iter_%d/PMs' % (last_idx+i+1), shell=True)
    subprocess.call('mv VPD*.png iter_%d' % (last_idx+i+1), shell=True)

    color_print('\tGenerando nuevo archivo de refstars', 'cyan')
    ids, pmx, pmy, nf, pmex, pmey = np.genfromtxt('iter_%d/PM_final.dat' % (last_idx+i+1), unpack=True, usecols=(0,3,4,6,8,9))
    pmr = np.sqrt(pmx**2 + pmy**2)
    pme = np.sqrt(np.square(pmex) + np.square(pmey))

    if os.path.isfile('zelimchisha'):
        rej_ids = np.genfromtxt('zelimchisha', unpack=True, usecols=(0,))
        id_mask = ~np.in1d(ids, rej_ids)
    else:
        id_mask = np.ones(len(ids)).astype(bool)

    mask = (pmr <= radio) * (nf >= nframes) * (pme <= max_err) * (id_mask)

    data = np.genfromtxt('iter_%d/PM_final.dat' % (last_idx+i+1))
    data = data[mask]

    fmt = '%d %.6f %.6f %.6f %.6f %.3f %d %d %.6f %.6f %.0f %.2f'
    hdr = 'ID RA DEC PM_X PM_Y MAG_K NFRAMES CFRAMES PMXE PMYE NEI NEI_STD'

    np.savetxt(refstars, data, fmt=fmt, header=hdr)

subprocess.call('rm -r PMs', shell=True)
subprocess.call('rm %s' % refstars, shell=True)

print '\nDone!'
