from __future__ import division, print_function
from astropy.utils.console import ProgressBar, color_print
import numpy as np
import pm_funcs
import argparse
import sys
import os

#Parametros y argumentos
parser = argparse.ArgumentParser(description='Tlineal 2.0')
parser.add_argument('INPUT', help='Archivo .mat')
#parser.add_argument('--max-mag', type=float, default=20.0, help='Corte superior en magnitud (Default 20)')
#parser.add_argument('--min-mag', type=float, default=8.0, help='Corte inferior en magnitud (Default 8)')
#parser.add_argument('--max-err', type=float, default=1.0, help='Maximo error a considerar (Default 1)')
#parser.add_argument('--comp', type=int, default=3, help='Nro de componentes para el Gaussian Mixture (Default 3)')

args   = parser.parse_args()
infile = args.INPUT
nepoca = int(infile.split('-')[1].split('_')[0])


nrefstars_tl, min_nei, rad_int, rad_ext, output, refer, sort_mag, \
local, ma1, ma2, mr1, mr2, mp1, mp2, rad_ref, x0, y0, lim, plot_ep, plot_del_ep, plot_del_xy, nprocs = pm_funcs.get_tlineal()
nrefstars = nrefstars_tl

color_print('TLINEAL: %s' % infile, 'lightcyan')
print('\tEpoca: %03d' % nepoca)
color_print('Parametros', 'yellow')
print('\tNro de refstars a considerar: %d' % nrefstars)
print('\tMin de vecinos a considerar:  %d' % min_nei)
print('\tRadio interior para buscar:   %d' % rad_int)
print('\tRadio exterior para buscar:   %d' % rad_ext)

#Chequea que exista archivo con refstars
if not os.path.isfile(refer):
    color_print('\nArchivo con estrellas de referencia no encontrado! --> %s' % refer, 'lightred')
    sys.exit(1)
