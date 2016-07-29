import os
import sys
import glob
import numpy as np
import pm_funcs
import argparse
import matplotlib.pyplot as plt

stilts_folder = os.path.dirname(os.path.realpath(__file__))

parser = argparse.ArgumentParser(description='CMD Bins')
parser.add_argument('K', type=int, help='Numero de epoca en Ks')
parser.add_argument('J', type=int, help='Numero de epoca en J')
parser.add_argument('--refstars', '-r', type=str, default=None, help='Destaca las estrellas en el archivo de refstars indicado')
parser.add_argument('--show-plot', '-p', action='store_true', help='Muestra el plot del CMD')

args = parser.parse_args()

cmd_modo, match, col1, col2, mag1, mag2, cmd_pdf = pm_funcs.get_CMD()

#Numero de epocas a usar
k_epoch = args.K
j_epoch = args.J

#Refstars
refst = args.refstars

if match=='RA DEC':
    match   = "RA DEC"
    omodo   = 'RD'

matcher = 'sky params=0.3'
if match == 'ID':
    matcher = 'exact'
    omodo   = 'ID'

output = 'CMD_%03d-%03d-%s.dat' % (k_epoch, j_epoch, omodo)

#Busca los catalogos
k_catalog = glob.glob('*k*%03d.datgc' % k_epoch)[0]
j_catalog = glob.glob('*j*%03d.datgc' % j_epoch)[0]

#Realiza el match usando K como referencia (1 and 2)
print 'Archivos para generar el CMD: %s %s\n' % (k_catalog, j_catalog)

os.system('java -jar %s/stilts.jar tmatch2 \
               in1=%s values1="%s" ifmt1=ascii \
               in2=%s values2="%s" ifmt2=ascii \
               matcher=%s find=best join=all1 \
               out=%s ofmt=ascii' % (stilts_folder,k_catalog,match,j_catalog,match,matcher,output) )

    #idk, rak, dek, xk, yk, magk, magek = np.genfromtxt(k_catalog, unpack=True)
    #idj, raj, dej, xj, yj, magj, magej = np.genfromtxt(j_catalog, unpack=True)

if cmd_pdf:
    ids, K, J = np.genfromtxt(output,unpack=True,usecols=(0,5,12))

    fig, ax = plt.subplots(figsize=[4,6])
    ax.plot(J-K,K,'.k',ms=1, alpha=1, rasterized=False)
    ax.set_xlabel('$J-K$')
    ax.set_ylabel('$K$')
    if cmd_modo=='manual':
        ax.set_xlim(col1,col2)
        ax.set_ylim(mag1,mag2)
    else:
        ax.invert_yaxis()
    if refst is not None:
        idr = np.genfromtxt(refst, unpack=True, usecols=(0,))
        idx = np.in1d(ids, idr)
        ax.plot((J-K)[idx], K[idx], '.r', ms=1, alpha=1, rasterized=False)

    print '\nGuardando PNG...'
    fig.savefig(output.replace('.dat', '.png'), dpi=200, bbox_inches='tight')
    if args.show_plot:
        plt.show()
print '\nDone!'
