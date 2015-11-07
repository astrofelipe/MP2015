#Crea un master frame usando como referencia(s) epocas dadas por el usuario
import os
import numpy as np
import glob
import sys
from astropy import wcs
from astropy.utils.console import ProgressBar, color_print
from astropy.io import fits
from scipy.optimize import curve_fit
from MP2015utils import *

#Para ignorar los warnings de WCS
import warnings
warnings.filterwarnings("ignore")

# PARAMETROS

iteraciones = 4     #Numero de iteraciones
threshold   = .5    #Porcentaje (0-1) en que debe aparecer una estrella para ser considerada
match_tol   = .3    #Tolerancia en arcsec para hacer los matches (0.3 arcsec ~ 1 pix en VVV)

# CARPETAS
radec_folder  = 'RADEC'  #Carpeta donde guarda los archivos con RA DEC
match_folder  = 'MATCH'  #Carpeta donde guarda los matches
stilts_folder = os.path.dirname(os.path.realpath(__file__)) #STILTS debe estar con el archivo .py


#####Informacion extra
if len(sys.argv) == 1:

    print
    ee = color_print('Como ejecutar','yellow')
    print 'python', sys.argv[0], 'path/to/catalogs/'
    print
    print 'Archivo "\033[1;33mzinfo_img\033[0m" debe estar en la carpeta de los datos'
    print 'Se requiere un archivo "\033[1;33mref_master\033[0m" que contenga la(s) imagen(es) de referencia'
    print 'Outputs:'

    sys.exit(1)

#Ruta de los .dao, .fits y zinfo
data_folder = sys.argv[1]
if not data_folder.endswith('/'):
    data_folder += '/'

#Crea las carpetas para guardar RADEC y matches
makedir(radec_folder)
makedir(match_folder)

#Lee ref_master y zinfo
color_print('Leyendo archivo con referencias...', 'cyan')
epn         = np.genfromtxt(data_folder+'zinfo_img', usecols=(0,), dtype='string')
se, el, yr  = np.genfromtxt(data_folder+'zinfo_img', unpack=True, usecols=(4,5,6))

km = np.array(['k' in e for e in epn])

refs = np.genfromtxt('ref_master', dtype='string')
kr   = refs[np.array(['k' in r for r in refs])]
jr   = refs[np.array(['j' in r for r in refs])]
print '\t%d referencias a usar' % len(kr)

#Busca archivos
color_print('Buscando archivos...', 'cyan')
k_files = np.sort(glob.glob('%s*k*.dao' % data_folder))
j_files = np.sort(glob.glob('%s*j*.dao' % data_folder))
print '\t%d Ks encontradas' % len(k_files)
print '\t%d J encontradas' % len(j_files)

#Convierte a RADEC
color_print('Convirtiendo XY a RADEC...','cyan')
print '\tConvirtiendo Ks'
ProgressBar.map(XYtoRADEC,k_files,multiprocess=True)
print '\tConvirtiendo J'
ProgressBar.map(XYtoRADEC,j_files,multiprocess=True)

#Crea el CMD
color_print('Creando CMDs...', 'cyan')
k_dats = np.sort(glob.glob('RADEC/*k*.dat'))
j_dats = np.sort(glob.glob('RADEC/*j*.dat'))

for k in k_dats:
    make_CMD(k, j_dats[0])
