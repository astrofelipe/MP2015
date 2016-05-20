from astropy.table import Table
import sys
import subprocess

#Input y output
inf  = sys.argv[1]
outf = sys.argv[1].replace('h5', 'dat')

#Abre el h5 como tabla
print '\nAbriendo archivo...'
t = Table.read(inf, path='data')
print '\nGuardando como ASCII (CSV)...'
subprocess.call("sed -i 's/nan//g' %s" % output, shell=True)
