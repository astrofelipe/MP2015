import os
import glob
import sys
from astropy.utils.console import ProgressBar

######### INPUT USUARIO
#python match1a1.py folder ref_cat
folder    = sys.argv[1] #path al directorio de las imagenes/catalogos
ref_cat   = folder+sys.argv[2] #catalogo de referencia
######### FIN INPUT USUARIO

catalog  = sorted(glob.glob(folder+'*k*.dat'))

#for cat in catalog:
def match(cat):

    if cat == ref_cat:
        continue

    #os.system('java -jar -Xmx4096M stilts.jar tmatch2 in1='+ref_cat+' values1="RA DEC" ifmt1=ascii in2='+cat+' values2="RA DEC" ifmt2=ascii matcher=sky params="0.3" find=best join=1and2 out='+cat.replace('.dat','_'+cat[-7:-4]+'.match')+' ofmt=ascii')
    os.system('java -jar stilts.jar tmatch2 in1='+ref_cat+' values1="RA DEC" ifmt1=ascii in2='+cat+' values2="RA DEC" ifmt2=ascii matcher=sky params="0.3" find=best join=1and2 out='+cat.replace('.dat','_'+cat[-7:-4]+'.match')+' ofmt=ascii')

ProgressBar.map(match,catalog,multiprocess=True)
