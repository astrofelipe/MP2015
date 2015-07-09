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

#Ejecutar como python tlineal_1a1.py <carpeta de las imagenes> <comodines de busqueda (ej *k*001.match)>
#Requisito: Archivo 'zinfo_img' en la carpeta de los datos

## PARAMETROS ##

vecinos = 51							#1 + vecinos (maximo)
rad_int  = 0							#Radio interior
rad_ext  = 0							#Radio exterior (Ambos en 0 = vecinos mas cercanos)
output  = 'test.pdf'					#PDF de Output (verificar siempre o se reemplazara sin avisar!)
refer   = '1116gc.dat'			#Catalogo con las estrellas locales
sort_mag = True							#Ordenar vecinos segun magnitud (toma los mas brillantes)
local   = True							#True para usar transf locales, False para usar global
ma1		= 11							#Limite inf de magnitudes para considerar los datos
ma2		= 15							#Limite sup de magnitudes para considerar los datos
ml1,ml2 = 11,15							#Corte en magnitud para estrellas locales

rad_clu = 300							#Radio en pixeles del cumulo
x0,y0 	= 1352,554						#Coordenadas del centro del cumulo en pixeles
mc1,mc2 = 11,14							#Limitde de magnitudes para plotear las estrellas del cumulo

lim 	= 8								#Limites del plot (cuadrado, por eso es uno)


#Con dist:  x0 1352 y0 554 r 280
#Libralato: x0 1795 y0 -6355 r 300 

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

def linear2(x,y,a,b,c):
	return a*x + b*y + c

delta_x = []
delta_y = []

def quadratic(coords,a,b,c,d,e,f):
	x,y = coords
	return a + b*x + c*y + d*np.square(x) + e*np.multiply(x,y) + f*np.square(y)

bid = np.genfromtxt(folder+refer,unpack=True,usecols=(0,))
print 'Locales: %d' % bid.size

fig_delta, ax_delta = plt.subplots(nro_rows*2,ncols=3,figsize=[5*3,2*nro_rows])
ad = np.ravel(ax_delta)

fig, ax = plt.subplots(nrows=nro_rows,ncols=3,figsize=[3.5*3,3.5*nro_rows])
fig.tight_layout()

for i,a in enumerate(np.ravel(ax)):
	if i == nro_arch:
		break

	print '\nIteracion %d/%d'%(i+1,ax.size)

	bid = np.genfromtxt(folder+refer,unpack=True,usecols=(0,))
	iid,x1,y1,m,id2,x2,y2,sep = np.genfromtxt(archivos[i],unpack=True,usecols=(0,3,4,5,7,10,11,14))
	#x1,y1,x2,y2 = 10000*np.array([x1,y1,x2,y2])
	mask			 = (m<ma2)*(m>ma1)
	#iid,x1,y1,m,x2,y2   = np.transpose([iid,x1,y1,m,x2,y2])[mask].T
	print 'Estrellas en Epoca: %d'%iid.size

	#print id[:-10]
	#print linear2(x2,y2,1,0,-17.24)[:-10]
	#print linear2(x2,y2,0,1,8.7682)[:-10]

	'''

	print linear2(2015.206,2070.265,0.999867313,-0.00146659,-17.2477)
	print linear2(2015.206,2070.265,0.001487624,0.999924547,8.7682)

	tx = linear2(x2,y2,0.999867313,-0.00146659,-17.2477)

	aaa = x1-tx
	ty = linear2(x2,y2,0.001487624,0.999924547,8.7682)
	eee = y1-linear2(x2,y2,0.001487624,0.999924547,8.7682)
	print aaa.min(),eee.min()
	print aaa.max(),eee.max()

	mask2 = aaa>3
	print mask2.sum()


	print id[mask2]
	print id2[mask2]
	print x1[mask2]
	print tx[mask2]
	print y1[mask2]
	print ty[mask2]

	fig, ax = plt.subplots(nrows=2)
	#ax.hist(aaa,bins=np.linspace(-4,4,100))
	#ax[0].scatter(x1,aaa,s=5,lw=0)
	#ax[1].scatter(y1,eee,s=5,lw=0)
	#ax.scatter(aaa,eee,lw=0,s=2,alpha=.5)
	#ax.set_aspect('equal')
	sca = ax[0].scatter(x1,y1,lw=.25,c=aaa,alpha=.5)
	scb = ax[1].scatter(x1,y1,lw=.25,c=eee,alpha=.5)
	plt.colorbar(sca,ax=ax[0])
	plt.colorbar(scb,ax=ax[1])
	plt.show()

	sys.exit()
	'''

	#Filtro por magnitud de las locales
	epinbu = np.in1d(iid,bid)
	buinep = np.in1d(bid,iid)

	bid,bx1,by1,bm,bx2,by2 = np.transpose([iid,x1,y1,m,x2,y2])[epinbu].T
	bid,bx1,by1,bm,bx2,by2 = np.transpose([bid,bx1,by1,bm,bx2,by2])[(bm>ml1)*(bm<ml2)].T

	epxy = np.transpose([bx2,by2])

	#Filtro por magnitud del catalogo
	iid,x1,y1,m,x2,y2 = np.transpose([iid,x1,y1,m,x2,y2])[mask].T

	print "Locales filtradas por mag: %d" % bid.size
	print "Estrellas del catalogo filtrado por mag: %d" %iid.size

	if local:

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
				guessx = [1,0,1]
				guessy = [0,1,1]

				poptx, pcovx = curve_fit(linear,coords,bx1[idx[k]],ftol=1e-10,xtol=1e-10)#,sigma=1e-3,absolute_sigma=True,p0=guessy)
				popty, pcovy = curve_fit(linear,coords,by1[idx[k]],ftol=1e-10,xtol=1e-10)#,sigma=1e-3,absolute_sigma=True,p0=guessx)
				
				pxt[k] += poptx
				pyt[k] += popty

				ctx[k] += linear([x2[k],y2[k]],*poptx)
				cty[k] += linear([x2[k],y2[k]],*popty)

				bar.update()

		figa,axa = plt.subplots(nrows=2,ncols=3)
		for i in range(3):
			axa[0,i].hist(pxt[:,i],bins=100)
			axa[1,i].hist(pyt[:,i],bins=100)
		plt.show()

	
	#Global
	else:
		poptx, pcovx = curve_fit(linear,[x2,y2],x1)
		popty, pcovy = curve_fit(linear,[x2,y2],y1)

		print poptx
		print popty

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
fig.savefig(output,dpi=200)
fig_delta.savefig('delta_'+output,dpi=200)

fig, ax = plt.subplots(nrows=2,figsize=[3*2,4*2])

for i in range(len(delta_x)):
	equis = np.zeros(len(delta_x[i])) + yrs[i]
	ax[0].scatter(equis,delta_x[i],c='#FF5500',lw=0,s=1,rasterized=True)
	ax[1].scatter(equis,delta_y[i],c='#0055FF',lw=0,s=1,rasterized=True)
fig.savefig('final_'+output,dpi=200)