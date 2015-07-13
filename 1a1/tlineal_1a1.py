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

#####Informacion extra
if len(sys.argv) == 1:

    print
    print 'Como ejecutar:', sys.argv[0]
    print 'python', sys.argv[0], 'path/to/catalogs/', '<texto para busqueda>'
    print
    print 'Archivo "zinfo_img" debe estar en la carpeta de los datos'
    print 'Output: archivo final_test.pdf de VPD para imagenes seleccionadas con respecto a archivo de referencia'
    print

    sys.exit(1)
#Ejecutar como python tlineal_1a1.py <carpeta de las imagenes> <comodines de busqueda (ej *k*001.match)>
#Requisito: Archivo 'zinfo_img' en la carpeta de los datos

## PARAMETROS ##

vecinos = 56							#1 + vecinos (maximo)
rad_int  = 1							#Radio interior
rad_ext  = 300							#Radio exterior (Ambos en 0 = vecinos mas cercanos)
output  = 'test.pdf'					#PDF de Output (verificar siempre o se reemplazara sin avisar!)
refer   = 'rgb.dat'			#Catalogo con las estrellas locales
sort_mag = True							#Ordenar vecinos segun magnitud (toma los mas brillantes)
local   = True							#True para usar transf locales, False para usar global
ma1, ma2 = 11,15							#Limites de magnitudes para considerar los datos
ml1,ml2 = 11,14							#Corte en magnitud para estrellas locales

rad_clu = 280							#Radio en pixeles del cumulo
x0,y0 	= 1352,554						#Coordenadas del centro del cumulo en pixeles
mc1,mc2 = 11,15							#Limitde de magnitudes para plotear las estrellas del cumulo

lim 	= 2								#Limites del plot (cuadrado, por eso es uno)

plot_PM = True						#Plot de delta x o delta y vs epocas

#Codigo

plt.style.use('ggplot')

folder   = sys.argv[1]

if sys.argv[2]=='-f':
	archivos = np.genfromtxt('files_match',unpack=True,dtype='string')

else:
	keywords = sys.argv[2]
	archivos = sorted(glob.glob(folder+keywords))

nro_arch = np.size(archivos)
nro_rows = nro_arch/3 + 1

se,el,yr = np.genfromtxt(folder+'zinfo_img',unpack=True,usecols=(4,5,6),skiprows=6)

def linear(coords,a,b,c):
	x,y = coords
	return a*x + b*y + c

def quadratic(coords,a,b,c,d,e,f):
	x,y = coords
	return a + b*x + c*y + d*np.square(x) + e*np.multiply(x,y) + f*np.square(y)

delta_x = []
delta_y = []

bid = np.genfromtxt(folder+refer,unpack=True,usecols=(0,))
print 'Locales: %d' % bid.size

fig_delta, ax_delta = plt.subplots(nro_rows*2,ncols=3,figsize=[5*3,2*nro_rows])
ad = np.ravel(ax_delta)

fig, ax = plt.subplots(nrows=nro_rows,ncols=3,figsize=[3.5*3,3.5*nro_rows])

for i,a in enumerate(np.ravel(ax)):
	if i == nro_arch:
		break

	print '\nIteracion %d/%d'%(i+1,ax.size)

	iid,x1,y1,m,id2,x2,y2,sep = np.genfromtxt(archivos[i],unpack=True,usecols=(0,3,4,5,7,10,11,14))

	print 'Estrellas en Epoca: %d'%iid.size

	#Filtro por magnitud de las locales
	bid    = np.genfromtxt(folder+refer,unpack=True,usecols=(0,))
	epinbu = np.in1d(iid,bid)

	bid,bx1,by1,bm,bx2,by2 = np.transpose([iid,x1,y1,m,x2,y2])[epinbu].T
	bid,bx1,by1,bm,bx2,by2 = np.transpose([bid,bx1,by1,bm,bx2,by2])[(bm>ml1)*(bm<ml2)].T

	epxy = np.transpose([bx2,by2])

	#Filtro por magnitud del catalogo
	mask = (m<ma2)*(m>ma1)
	iid,x1,y1,m,x2,y2 = np.transpose([iid,x1,y1,m,x2,y2])[mask].T

	print "Locales filtradas por mag: %d" % bid.size
	print "Estrellas del catalogo filtrado por mag: %d" %iid.size

	if local:

		nbrs = NN(n_neighbors=vecinos, algorithm='auto').fit(epxy)

		if rad_ext!=0:
			dist, idx = nbrs.radius_neighbors(np.transpose([x2,y2]),radius=rad_ext)
			nbors 	  = np.array([len(d) for d in dist])

			mednbors, minnbors = np.median(nbors),nbors.min()

			for j in range(len(dist)):
				#Elimina la misma estrella y las interiores a rad_int
				msk	    = (dist[j]!=0)*(dist[j]>rad_int)
				dist[j] = dist[j][msk]
				idx[j]  = idx[j][msk]

				#Toma las mas brillantes
				if len(dist[j])>vecinos:
					if sort_mag:
						lbm		= bm[idx[j]]
						midx 	= np.argsort(lbm)[:vecinos]
						dist[j] = dist[j][midx]
						idx[j]  = idx[j][midx]

				#Toma las mas cercanas
					else:
						midx 	= np.argsort(dist[j])[:vecinos]
						dist[j] = dist[j][midx]
						idx[j]  = idx[j][midx]

		#Vecinos mas cercanos
		else:
			dist,idx = nbrs.kneighbors(np.transpose([x2,y2]))
			idx		 = idx[:,1:]
			dist 	 = dist[:,1:]

		means = np.array([np.mean(d) for d in dist])
		nbors 	  = np.array([len(d) for d in dist])
		mednbors, minnbors = np.median(nbors),nbors.min()

		ctx = np.zeros(x1.size)
		cty = np.zeros(y1.size)

		pxt = np.zeros((x1.size,3))
		pyt = np.zeros((y1.size,3))

		with ProgressBar(x1.size) as bar:
			for k in range(x1.size):
				coords = np.transpose([bx2,by2])[idx[k]].T
				#coords = (coords.T - np.mean(coords,axis=1)).T

				guessx = [1,0,1]
				guessy = [0,1,1]

				ep1_x = bx1[idx[k]]# - np.mean(bx1[idx[k]])
				ep1_y = by1[idx[k]]# - np.mean(by1[idx[k]])

				poptx, pcovx = curve_fit(linear,coords,ep1_x)
				popty, pcovy = curve_fit(linear,coords,ep1_y)

				pxt[k] += poptx
				pyt[k] += popty

				ctx[k] += linear([x2[k],y2[k]],*poptx)
				cty[k] += linear([x2[k],y2[k]],*popty)

				bar.update()

	#Global
	else:
		poptx, pcovx = curve_fit(linear,[x2,y2],x1)
		popty, pcovy = curve_fit(linear,[x2,y2],y1)

		ctx = linear([x2,y2],*poptx)
		cty = linear([x2,y2],*popty)

		means = -1
		nbors = -1
		mednbors, minnbors = -1,-1


	#MAIN PLOT

	m_clu = (m>mc1)*(m<mc2)
	clust = (np.sqrt((x1-x0)**2 + (y1-y0)**2) < rad_clu)

	a.scatter((x1-ctx)[~clust],(y1-cty)[~clust],s=1,rasterized=True,edgecolor='',color='#0055FF',lw=.5)
	a.scatter((x1-ctx)[clust*m_clu],(y1-cty)[clust*m_clu],s=1.25,rasterized=True,edgecolor='',color='#FF5500',lw=.5)
	a.text(.1,.1,u'$S = %f$\n$E = %s$\n$N = %d/%d$\n$Med_d, Max_d = %.3f, %.3f$\n$Med_n,Min_n = %d,%d$'%(se[i+1],el[i+1],(clust*m_clu).sum(),clust.size,np.median(means),np.max(means),mednbors, minnbors),transform = a.transAxes,alpha=.66,fontsize=10)
	a.text(.9,.9,u'$%d$'%(i+2),transform = a.transAxes,alpha=.66,fontsize=14)

	xmean_b = np.mean((x1-ctx)[~clust])
	ymean_b = np.mean((y1-cty)[~clust])
	xstd_b	= np.std((x1-ctx)[~clust])
	ystd_b	= np.std((y1-cty)[~clust])

	xmean_r = np.mean((x1-ctx)[clust*m_clu])
	ymean_r = np.mean((y1-cty)[clust*m_clu])
	xstd_r	= np.std((x1-ctx)[clust*m_clu])
	ystd_r	= np.std((y1-cty)[clust*m_clu])

	a.text(.1,.83,u'M = $%.4f / %.4f$\nS = $%.3f / %.3f$' % (xmean_r,ymean_r,xstd_r,ystd_r),color='#FF5500',transform=a.transAxes)
	a.text(.1,.69,u'M = $%.4f / %.4f$\nS = $%.3f / %.3f$' % (xmean_b,ymean_b,xstd_b,ystd_b),color='#0055FF',transform=a.transAxes)
	a.set_xlim(-lim,lim)
	a.set_ylim(-lim,lim)
	a.set_aspect('equal')

	#x vs dx / y vs dy

	dx = x1-ctx
	dy = y1-cty

	print '        dx, dy'
	print 'max:    %f, %f' % (dx.max(),dy.max())
	print 'mean:   %f, %f' % (np.mean(dx),np.mean(dy))
	print 'median: %f, %f' % (np.median(dx),np.median(dy))
	print 'std:    %f, %f' % (np.std(dx),np.std(dy))

	ad[2*i].scatter(x1[~clust],dx[~clust],s=1,rasterized=True,lw=0,color='#0055FF')
	ad[2*i].scatter(x1[clust*m_clu],dx[clust*m_clu],s=1,rasterized=True,lw=0,color='#FF5500')

	ad[2*i+1].scatter(y1[~clust],dy[~clust],s=1,rasterized=True,lw=0,color='#0055FF')
	ad[2*i+1].scatter(y1[clust*m_clu],dy[clust*m_clu],s=1,rasterized=True,lw=0,color='#FF5500')

	#"FINAL" PLOT
	delta_x.append(dx)
	delta_y.append(dy)

	if (i%5)==0:
		print '\nGuardando...\n'
		fig.savefig(output,dpi=200)
		fig_delta.savefig('delta_'+output,dpi=200)
fig.tight_layout()
fig.savefig(output,dpi=200)
#fig_delta.savefig('delta_'+output,dpi=200)

#Final Plot cont
if plot_PM:
	yrs = (yr-yr[0])/365.25

	fit_meansx = np.array([np.mean(dd) for dd in delta_x])
	fit_meansy = np.array([np.mean(dd) for dd in delta_y])

	fig, ax = plt.subplots(nrows=2,figsize=[3*2,4*2])

	for i in range(len(delta_x)):
		equis = np.zeros(len(delta_x[i])) + yrs[i]
		ax[0].scatter(yrs[:len(delta_x)],fit_meansx,c='#FF5500',lw=0,s=1,rasterized=True)
		ax[1].scatter(yrs[:len(delta_y)],fit_meansy,c='#0055FF',lw=0,s=1,rasterized=True)

	coeffx = np.polyfit(yrs[:len(delta_x)],fit_meansx,1)
	coeffy = np.polyfit(yrs[:len(delta_y)],fit_meansy,1)

	ax[0].plot(yrs[:len(delta_x)],np.polyval(coeffx,yrs[:len(delta_x)]))
	ax[1].plot(yrs[:len(delta_y)],np.polyval(coeffy,yrs[:len(delta_y)]))

	fig.savefig('final_'+output,dpi=200)
