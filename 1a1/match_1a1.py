import os
import glob
import sys
import multiprocessing as mp
import numpy as np

#####Informacion extra
if len(sys.argv) == 1:

    print
    print 'Como ejecutar:', sys.argv[0]
    print 'python', sys.argv[0], 'input list', 'reference_catalog'
    print
    print 'Output: archivos *.mat'
    print 'Match con el catalogo de referencia de cara archivo *k* en el path'
    print

    sys.exit(1)

#PARAMETROS
from zparams_pm import modo, tol

######### INPUT USUARIO
#python match1a1.py folder ref_cat
#folder    = sys.argv[1] #path al directorio de las imagenes/catalogos
inputs    = sys.argv[1]
ref_cat   = sys.argv[2] #catalogo de referencia
######### FIN INPUT USUARIO

print 'Iniciando con...'
print 'modo: %s' % modo
print 'tol:  %.f' % tol

stilts_folder = os.path.dirname(os.path.realpath(__file__))
catalog  = np.genfromtxt(inputs, dtype='string')
#catalog  = sorted(glob.glob(folder+'*k*.dat'))

matcher = 'sky params=%.3f' % tol
if modo=='ID':
    matcher = 'exact'

#for cat in catalog:
def match(cat):

    if cat == ref_cat:
        return
    zeros = ref_cat.split('.')[0][-3:]
    #os.system('java -jar -Xmx4096M stilts.jar tmatch2 in1='+ref_cat+' values1="RA DEC" ifmt1=ascii in2='+cat+' values2="RA DEC" ifmt2=ascii matcher=sky params="0.3" find=best join=1and2 out='+cat.replace('.dat','_'+cat[-7:-4]+'.match')+' ofmt=ascii')
    execute = 'java -jar '+ stilts_folder+'/stilts.jar tmatch2 in1='+ref_cat+' values1='+modo+' ifmt1=ascii in2='+cat+' values2='+modo+' ifmt2=ascii matcher='+matcher+' find=best join=1and2 out='+cat.replace('.dat','_'+zeros+'.mat')+' ofmt=ascii'
    #print execute,'\n'
    os.system(execute)

cpus = (1 + mp.cpu_count()/2)

pool = mp.Pool(processes=cpus)
pool.map(match,catalog)

#ProgressBar.map(match,catalog,multiprocess=True)
