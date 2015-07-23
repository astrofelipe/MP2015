import os
import sys
import glob
import numpy as np
from sklearn.neighbors import NearestNeighbors as NN
from astropy.utils.console import ProgressBar, color_print

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

color_print('Leyendo informacion de epocas')
se,el,yr = np.genfromtxt(folder+'zinfo_img',unpack=True,usecols=(4,5,6),skiprows=6)

color_print('Recopilando archivos de epocas...','cyan')
epochs = glob.glob('./%s/*.*' % matched_epochs)

color_print('Realizando match de la MF con las epocas','cyan')
ejecuta = 'java -jar %s/stilts.jar tmatch2 in1=%s values1="ID" ifmt1=ascii ' % (stilts_folder, master)

for i in range(len(epochs)):
def mf_match(ep):
    ej2 = 'in2=%s values1="ID" ifmt2=ascii matcher=exact find=best join=1and2 out=./%s/%s ofmt=ascii progress=none' % (ep, match_master, ep.replace('.match','.mfma'))

ProgressBar.map(mf_match,epochs,multiprocess=True)
