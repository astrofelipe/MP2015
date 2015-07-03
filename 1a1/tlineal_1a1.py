import numpy as np
import os
import sys
import glob
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit
from scipy.spatial import cKDTree
from sklearn.neighbors import NearestNeighbors as NN 
from sklearn.neighbors import RadiusNeighborsClassifier as RN 
from astropy.utils.console import ProgressBar

#Ejecutar como python tlineal_1a1.py <carpeta con los *_ddd.match> con ddd nro de epoca de referencia
#Requisito: Archivo 'zinfo_img' en la carpeta de los datos

## PARAMETROS ##

vecinos = 11							#1 + vecinos (maximo)
rad_int  = 0							#Radio interior
rad_ext  = 500							#Radio exterior (Ambos en 0 = vecinos mas cercanos)
sort_mag = True							#Ordenar vecinos segun magnitud (toma los mas brillantes)
output  = 'test.pdf'					#PDF de Output (verificar siempre o se reemplazara sin avisar!)
refer   = '1114.dat'			#Catalogo con las estrellas locales
ma1		= 11							#Limite inf de magnitudes para considerar los datos
ma2		= 15							#Limite sup de magnitudes para considerar los datos

rad_clu = 280							#Radio en pixeles del cumulo
x0,y0 	= 1352,554 						#Coordenadas del centro del cumulo en pixeles
mc1,mc2 = 11,14							#Limitde de magnitudes para plotear las estrellas del cumulo

#Codigo

plt.style.use('ggplot')

folder   = sys.argv[1]
keywords = sys.argv[2]
archivos = sorted(glob.glob(folder+keywords))
nro_arch = np.size(archivos)
nro_rows = nro_arch/3 + 1

se,el,yr = np.genfromtxt(folder+'zinfo_img',unpack=True,usecols=(4,5,6),skiprows=6)
yrs = (yr-yr[0])/365.25

def linear(coords,a,b,c):
	x,y = coords
	return a*x + b*y + c

def quadratic(coords,a,b,c,d,e,f):
	x,y = coords
	return a + b*x + c*y + d*np.square(x) + e*np.multiply(x,y) + f*np.square(y)

bid,bx,by = np.genfromtxt(folder+refer,unpack=True,usecols=(0,3,4))

fig, ax = plt.subplots(nrows=nro_rows,ncols=3,figsize=[3.5*3,3.5*nro_rows])
fig.tight_layout()

for i,a in enumerate(np.ravel(ax)):
	if i == nro_arch:
		break

	print 'Iteracion %d/%d'%(i+1,ax.size)

	#r1,d1,m,r2,d2 = np.genfromtxt('%s.dat'%(70*j/8+(i+1)),unpack=True,usecols=(1,2,5,7,8))
	id,x1,y1,m,x2,y2 = np.genfromtxt(archivos[i],unpack=True,usecols=(0,3,4,5,10,11))
	mask			 = (m<ma2)*(m>ma1)
	id,x1,y1,m,x2,y2   = np.transpose([id,x1,y1,m,x2,y2])[mask].T

	epinbu = np.in1d(id,bid)
	buinep = np.in1d(bid,id)

	print id[epinbu].size,bid[buinep].size

	epxy = np.transpose([x2,y2])[epinbu]
	epm  = m[epinbu]

	nbrs = NN(n_neighbors=vecinos, algorithm='auto').fit(epxy)

	if rad_int!=0 and rad_ext!=0:
		dist, idx = nbrs.radius_neighbors(np.transpose([x2,y2]),radius=rad_ext)
		nbors 	  = np.array([len(d) for d in dist])

		mednbors, minnbors = np.median(nbors),nbors.min()

		for j in range(len(dist)):
			#Elimina la misma estrella
			msk	    = (dist[j]!=0)*(dist[j]>=rad_int)
			dist[j] = dist[j][msk]
			idx[j]  = idx[j][msk]


			#Toma las mas brillantes
			if len(dist[j])>vecinos:
				if sort_mag:
					midx 	= np.argsort(epm[idx[j]])[:vecinos]
					dist[j] = dist[j][midx]
					idx[j]  = idx[j][midx]
			
			#Toma las mas cercanas
				else:
					midx 	= np.argsort(dist[j])[:vecinos]
					dist[j] = dist[j][midx]
					idx[j]  = idx[j][midx]
			
	else:
		dist,idx = nbrs.kneighbors(np.transpose([x2,y2]))
		idx		 = idx[:,1:]
		dist 	 = dist[:,1:]

	means = np.array([np.mean(d) for d in dist])
	nbors 	  = np.array([len(d) for d in dist])
	mednbors, minnbors = np.median(nbors),nbors.min()

	ctx = np.zeros(x1.size)
	cty = np.zeros(y1.size)
	
	with ProgressBar(x1.size) as bar:
		for k in range(x1.size):
			coords = np.transpose([x2,y2])[idx[k]].T
			poptx, pcovx = curve_fit(linear,coords,x1[idx[k]])
			popty, pcovy = curve_fit(linear,coords,y1[idx[k]])
			
			the_x = x2[k]
			the_y = y2[k]

			ctx[k] += linear([the_x,the_y],*poptx)
			cty[k] += linear([the_x,the_y],*popty)

			bar.update()

	m_clu = (m>mc1)*(m<mc2)
	clust = (np.sqrt((x1-x0)**2 + (y1-y0)**2) < rad_clu)

	a.scatter((x1-ctx)[~clust],(y1-cty)[~clust],s=1,rasterized=True,edgecolor='',color='#0055FF',lw=.5)
	a.scatter((x1-ctx)[clust*m_clu],(y1-cty)[clust*m_clu],s=1.25,rasterized=True,edgecolor='',color='#FF5500',lw=.5)
	a.text(.1,.1,u'$S = %f$\n$E = %s$\n$N = %d/%d$\n$Med_d, Max_d = %.3f, %.3f$\n$Med_n,Min_n = %d,%d$'%(se[i+1],el[i+1],(clust*m_clu).sum(),clust.size,np.median(means),np.max(means),mednbors, minnbors),transform = a.transAxes,alpha=.66,fontsize=10)
	a.text(.1,.9,u'$%d$'%(i+2),transform = a.transAxes,alpha=.66,fontsize=14)
	#a.set_xlim(-6e-5,6e-5)
	#a.set_ylim(-6e-5,6e-5)
	a.set_xlim(-8,8)
	a.set_ylim(-8,8)
	a.set_aspect('equal')

	dx = x1-ctx
	dy = y1-cty

	print dx.max(),dy.max()
	print np.mean(dx),np.mean(dy)
	print np.median(dx),np.median(dy)

	#np.savetxt()

	#a.set_xticklabels([])
	#a.set_yticklabels([])

		#plt.show()
	if (i%5)==0:
		print 'Guardando...\n'
		#fig.suptitle('Neigh = %d'%vecinos)
		fig.savefig(output,dpi=200)
fig.savefig(output,dpi=200)