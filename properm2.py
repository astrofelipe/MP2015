import os
import sys
import glob
import numpy as np
from sklearn.neighbors import NearestNeighbors as NN
from astropy.utils.console import ProgressBar, color_print
from scipy.optimize import curve_fit

#PARAMETROS
vecinos = 56
output  = 'PM_.dat'
locales = 'disk.dat'
master  = 'Master.dat'
ma1,ma2 = 9, 20
ml1,ml2 = 11, 14

vc_pix = 0.34

match_folder = 'matched_epochs'
match_master = 'matched_mf'
pm_folder    = 'epochs_pm'

stilts_folder = os.path.dirname(os.path.realpath(__file__)) #STILTS debe estar con el archivo .py

#FUNCIONES
def makedir(directory):
	if not os.path.exists(directory):
		os.makedirs(directory)

def linear(coords,a,b,c):
	x,y = coords
	return a + b*x + c*y

#PIPELINE
folder = sys.argv[1]
makedir(match_folder)
makedir(match_master)

color_print('Leyendo informacion de epocas')
se,el,yr = np.genfromtxt(folder+'zinfo_img',unpack=True,usecols=(4,5,6),skip_header=6)

color_print('Recopilando archivos de epocas...','cyan')
epochs = glob.glob('./%s/*.*' % match_folder)

color_print('Realizando match de la MF con las epocas','cyan')
ejecuta = 'java -jar %s/stilts.jar tmatch2 in1=./%s values1="ID" ifmt1=ascii ' % (stilts_folder, master)

def mf_match(ep):
    ej2 = 'in2=%s values2="ID_1" ifmt2=ascii icmd2=\'keepcols "ID_1 X Y"\' matcher=exact find=best join=1and2 out=./%s/%s ofmt=ascii progress=none ocmd="delcols ID_1"' % (ep, match_master, ep.split('/')[-1].replace('.match','.mfma'))
    os.system(ejecuta + ej2)

ProgressBar.map(mf_match,epochs,multiprocess=True)

color_print('Realizando transformaciones lineales')
matches = glob.glob('./%s/*.*' % match_master)
bid     = np.genfromtxt(locales,unpack=True,usecols=(0,))

def shift(ep):
    ids,x1,y1,mk,mj,x2,y2 = np.genfromtxt(ep,unpack=True)

    local_mask = np.in1d(ids,bid)
    lid,lx1,ly1,lx2,ly2  = np.transpose([ids,x1,y1,x2,y2])[local_mask].T

    loc_xy = np.transpose([lx2,ly2])
    nbrs   = NN(n_neighbors=vecinos, algorithm='auto').fit(loc_xy)

    coords = np.transpose([x2,y2])
    dist, idx = nbrs.kneighbors(coords)
    idx  = idx[:,1:]
    dist = dist[:,1:]

    ctx = np.zeros(x1.size)
    cty = np.zeros(y1.size)

    for i in range(x1.size):
        star_locales = np.transpose(coords)[idx[i]].T
        ep1_x = bx1[idx[i]]
        ep1_y = by1[idx[i]]

        poptx, pcovx = curve_fit(linear,star_locales,ep1_x)
        popty, pcovy = curve_fit(linear,star_locales,ep1_y)

        ctx[i] += linear([x2[i],y2[i]],*poptx)
        cty[i] += linear([x2[i],y2[i]],*popty)

    shift_x = x1 - ctx
    shift_y = y1 - cty

    hdr  = 'ID X Y MAG_K MAG_J PMX PMY'
    fmt  = '%d %.3f %.3f %.3f %.3f %f %f'
    data = np.transpose([ids,x1,y1,mk,mj,shift_x,shift_y])
    np.savetxt(ep.replace('.mfma','.pm'),data,header=hdr,fmt=fmt)

ProgressBar.map(shift,matches,multiprocess=True)
