import numpy as np
from astropy import wcs
from astropy.io import fits
import os

def XYtoRADEC(ep):
    cfn = ep
    ffn = ep.replace('.dao','.fits')

    rid,rx,ry,rmag = np.genfromtxt(cfn,usecols=range(4),skip_header=3,unpack=True)

    hdr    = fits.getheader(ffn)
    w      = wcs.WCS(hdr)
    rcoo   = np.transpose([rx,ry])
    ra,dec = np.transpose(w.wcs_pix2world(rcoo,1))

    head = 'ID RA DEC X Y MAG'
    fmt  = '%d %f %f %.3f %.3f %.3f'
    of   = ep.split('/')[-1].replace('.dao','.dat')
    np.savetxt('./%s/%s' % (radec_folder,of),np.transpose([rid,ra,dec,rx,ry,rmag]),header=head,fmt=fmt)

def makedir(directory):
	if not os.path.exists(directory):
		os.makedirs(directory)

#Transformacion lineal
def linear(coords,a,b,c):
	x,y = coords
	return a + b*x + c*y
