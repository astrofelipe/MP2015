import numpy as np
from astropy import wcs
from astropy.io import fits as pf
from astropy import wcs
from astropy.utils.console import ProgressBar, color_print
import os
import sys

path = sys.argv[1]
os.chdir(path)

#Transforma de XY a RADEC
def XYtoRADEC(ep):
	ffn,cfn = ep

	id,x,y,mag = np.loadtxt(cfn,usecols=range(4),skiprows=3,unpack=True)
	id = id.astype(int)
	
	hdr    = pf.open(ffn)[0].header
	w      = wcs.WCS(hdr)
	ra,dec = np.transpose(w.wcs_pix2world(np.transpose([x,y]),1))
	
	head = 'ID RA DEC X Y MAG'
	fmt  = '%d %.7f %.7f %.3f %.3f %.3f'
	np.savetxt(ffn.replace('fits','dat'),np.transpose([id,ra,dec,x,y,mag]),header=head,fmt=fmt)

archivos  = np.sort([f for f in os.listdir('./')])
fits	  = np.sort([f for f in archivos if f.endswith('.fits')])
catalog	  = np.sort([f for f in archivos if f.endswith('.dao')])

color_print('-Obteniendo RADEC...','cyan')
ProgressBar.map(XYtoRADEC,np.transpose([fits,catalog]),multiprocess=True)

