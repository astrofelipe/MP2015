import numpy as np
from astropy import wcs
from astropy.io import fits as pf
from astropy import wcs
from astropy.utils.console import ProgressBar, color_print
import os
import sys
import warnings

warnings.filterwarnings("ignore")

#####Informacion extra
if len(sys.argv) == 1:

	print
	print 'Como ejecutar:', sys.argv[0]
	print 'python', sys.argv[0], 'path/to/catalogs/ <extension>'
	print '(extension sin punto)'
	print
	print 'Output: archivos *.dat con XY convertidos a RADEC sobreescribiendolos'
	print

	sys.exit(1)

path = sys.argv[1]
ext  = sys.argv[2]

#Transforma de XY a RADEC
def XYtoRADEC(ep):
	ffn,cfn = ep

	id,x,y,mag,err = np.loadtxt(path+cfn,usecols=range(5),skiprows=3,unpack=True)
	id = id.astype(int)

	hdr    = pf.open(path+ffn)[0].header
	w      = wcs.WCS(hdr)
	ra,dec = np.transpose(w.wcs_pix2world(np.transpose([x,y]),1))

	head = 'ID RA DEC X Y MAG ERR'
	fmt  = '%d %.7f %.7f %.3f %.3f %.3f %.3f'
	np.savetxt(ffn.replace('fits','dat'),np.transpose([id,ra,dec,x,y,mag,err]),header=head,fmt=fmt)

archivos  = np.sort([f for f in os.listdir(path)])
fits	  = np.sort([f for f in archivos if f.endswith('.fits')])
catalog	  = np.sort([f for f in archivos if f.endswith('.'+ext)])

color_print('-Obteniendo RADEC...','cyan')
ProgressBar.map(XYtoRADEC,np.transpose([fits,catalog]),multiprocess=True)
