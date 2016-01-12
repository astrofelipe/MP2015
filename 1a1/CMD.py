import os
import sys
import glob
import numpy as np
import pm_params

stilts_folder = os.path.dirname(os.path.realpath(__file__))

#####Informacion extra
if len(sys.argv) == 1:

    print
    print 'Como ejecutar:', sys.argv[0]
    print 'python', sys.argv[0], 'Numero de epoca en Ks', 'Numero de epoca en J'
    print
    print 'Output: CMD_*-*.dat match de las epocas en J y Ks '
    print 'Opcionalmente hace el plot del CMD si esta la opcion "-p" luego de la epoca en J'
    print

    sys.exit(1)

cmd_modo, col1, col2, mag1, mag2 = pm_params.get_CMD()

#Numero de epocas a usar
k_epoch = int(sys.argv[1])
j_epoch = int(sys.argv[2])

output = 'CMD_%03d-%03d.dat' % (k_epoch, j_epoch)

#Busca los catalogos
k_catalog = glob.glob('*k*%03d.datgc' % k_epoch)[0]
j_catalog = glob.glob('*j*%03d.datgc' % j_epoch)[0]

#Realiza el match usando K como referencia (1 and 2)
print 'Archivos para generar el CMD: %s %s\n' % (k_catalog, j_catalog)
os.system('java -jar %s/stilts.jar tmatch2 \
           in1=%s values1="RA DEC" ifmt1=ascii \
           in2=%s values2="RA DEC" ifmt2=ascii \
           matcher=sky params=0.3 find=best join=all1 \
           out=%s ofmt=ascii' % (stilts_folder,k_catalog,j_catalog,output) )

if len(sys.argv) >3:
    if sys.argv[3] == '-p':
        import matplotlib.pyplot as plt
        K,J = np.genfromtxt(output,unpack=True,usecols=(5,12))

        fig, ax = plt.subplots()
        ax.plot(J-K,K,'.k',ms=2,alpha=.5)
        ax.set_xlabel('$J-K$')
        ax.set_ylabel('$K$')
        if cmd_modo=='manual':
            ax.set_xlim(col1,col2)
            ax.set_ylim(mag1,mag2)
        else:
            ax.invert_yaxis()
        plt.show()
