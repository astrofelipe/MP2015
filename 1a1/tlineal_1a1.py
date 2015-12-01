import numpy as np
import os
import sys
import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.spatial import cKDTree
from sklearn.neighbors import NearestNeighbors as NN
from sklearn.neighbors import RadiusNeighborsClassifier as RN
from astropy.io import fits
from astropy.utils.console import ProgressBar
from astropy.stats.funcs import mad_std

#####Informacion extra
if len(sys.argv) == 1:

    print
    print 'Como ejecutar:', sys.argv[0]
    print 'python', sys.argv[0], 'path/to/catalogs/', '<texto para busqueda>'
    print
    print 'Archivo "zinfo_img" debe estar en la carpeta de los datos'
    print 'Funciona con una lista de catalogos (archivo files_match)'
    print 'si se cambia "<texto para busquedas>" por "-f"'
    print 'Outputs:'
    print 'output.psf: delta_x vs delta_y para cada catalogo'
    print 'output_del_xy.pdf: delta_x vs x + delta_y vs y'
    print 'output_del_ep.pdf: med_delta vs time'
    print 'El programa se interrumpira si encuentra menos de 4 estrellas locales.'

    sys.exit(1)

#Ejecutar como python tlineal_1a1.py <carpeta de las imagenes> <comodines de busqueda (ej *k*001.match)>
#Requisito: Archivo 'zinfo_img' en la carpeta de los datos

## PARAMETROS ##
## Las estrellas de referencia (refstars) son las que se usan para llevar a cabo las transformaciones
## de coordenadas de cada estrella analizada

nrefstars = 51         #numero de refstars(+1) deseadas
rad_int = 1         #Radio interior
rad_ext = 300         #Radio exterior (0 -> selecciona nrefstars mas cercanas)
output  = 'test1'     #PDF de Output (reemplazara sin avisar!)
refer   = 'refstars1.gc' #Catalogo con las estrellas de referencia
sort_mag = False      #Sort refstars segun mag y toma las nrefstars mas brillantes
local   = True           #True para usar transf locales, False para usar global
ma1,ma2 = 11.0,14.0      #Corte en magnitud para considerar estrellas a analizar
mr1,mr2 = 11.0,14.0      #Corte en magnitud para las refstars
mp1,mp2 = 11.0,12.0     #Corte en magnitud para plotear las refstars

rad_ref = 99999999       #Radio (pix) dentro del cual una refstar se considera para plot
x0,y0    = 1352,554       #Coordenadas centrales del circulo a considerar

lim    = 2         #Limites del plot (cuadrado, por eso es uno)

plot_del_ep = True     #Plot delta vs epocas
plot_del_xy = True      #Plot delta vs coor_x o coor_y

############################################# CODIGO

folder   = sys.argv[1]
if not os.path.isfile(folder+refer):
    print '\nNo encuentro archivo con estrellas de referencia -->', refer
    print 'Bye Bye...'
    sys.exit(1)

#if not os.path.isfile(folder+zinfo_img):
#    print '\nNo encuentro archivo con la info -->', zinfo_img
#    print 'Bye Bye...'
#    sys.exit(1)

if sys.argv[2]=='-f':
    #numpy usa genfromtxt para abrir archivos
    #numpy lee numeros, por eso en este caso se especifica string.
    archivos = np.genfromtxt('files_match',unpack=True,dtype='string')

else:
    keywords = sys.argv[2]
    #glog es un paquete, el segundo glob es como un ls
    archivos = sorted(glob.glob(folder+keywords))

nro_arch = np.size(archivos)
print 'archivos =',nro_arch
##esto es para asignar numero de filas al disegno del plot output.pdf
nro_rows  = nro_arch/3 + 1
nro_rows2 = nro_arch/2 + 1
print 'rows =',nro_rows
#?que hace este?
nro_epoca = np.sort([int(f.split('-')[1].split('_')[0]) for f in archivos])
print 'epocas =',nro_epoca

se,el,yr = np.genfromtxt(folder+'zinfo_img',unpack=True,usecols=(4,5,6))
zn = np.genfromtxt(folder+'zinfo_img', unpack=True, usecols=(0,), dtype='string')
ky = np.array(['k' in z for z in zn])
yr = yr[ky]

se = se[ky][nro_epoca-1]
el = el[ky][nro_epoca-1]

def linear(coords,a,b,c):
    x,y = coords
    return a*x + b*y + c

def quadratic(coords,a,b,c,d,e,f):
    x,y = coords
    return a + b*x + c*y + d*np.square(x) + e*np.multiply(x,y) + f*np.square(y)

def recta(x,a,b):
    return a*x + b

def makedir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

#Promedios de los 2 grupos (campo y refstars)
#Creo vectores con tantos ceros como numero de archivos. Los vectores son horiz.
#?por que hay que crearlos?
meansx_ref = np.zeros(nro_arch)
meansy_ref = np.zeros(nro_arch)
stdx_ref   = np.zeros(nro_arch)
stdy_ref   = np.zeros(nro_arch)

meansx_fie = np.zeros(nro_arch)
meansy_fie = np.zeros(nro_arch)
stdx_fie   = np.zeros(nro_arch)
stdy_fie   = np.zeros(nro_arch)
#print 'meansx_ref',meansx_ref
#sys.exit()

#fig_delta, ax_delta = plt.subplots(nro_rows*2,ncols=3,figsize=[5*3,2*nro_rows])
#con el que sigue plotea test_del_xy.pdf de 2 en 2
fig_delta, ax_delta = plt.subplots(nrows=nro_arch,ncols=2,figsize=[5*1.5,3*nro_rows])
#?que hace ravel?
ad = np.ravel(ax_delta)
#print 'ad', ad.shape
#sys.exit()
fig, ax = plt.subplots(nrows=nro_rows,ncols=3,figsize=[3.5*3,3.5*nro_rows])

#aca empieza el for para todos los catalogos de input
for i,a in enumerate(np.ravel(ax)):
    if i == nro_arch:
            break

    #?que es ax?
    print '\nIteracion %d/%d'%(i+1,ax.size)

    #abro el catalogo i
    id1,x1,y1,m,id2,x2,y2,sep = np.genfromtxt(archivos[i],unpack=True,usecols=(0,3,4,5,7,10,11,14))
    print '\nAnalizando el catalogo %s' %(archivos[i])
    print 'Estrellas en Epoca: %d' %id1.size

    #abro archivo de refstars para contar
    rid0 = np.genfromtxt(folder+refer,unpack=True,usecols=(0,))
    print 'Estrellas de ref en archivo refstars: %d' %(rid0.size)
    #np.in1d es un match en 1 dimension y genera un vector de tantos true or false
    #como entradas hay en catalogo i (manda el primero)
    epinbu = np.in1d(id1,rid0)
    #print 'epinbu', epinbu.size
    #sys.exit()

    #selecciono solo las estrellas del catalogo en cuestion que son refstars
    rid1,rx1,ry1,rm,rx2,ry2 = np.transpose([id1,x1,y1,m,x2,y2])[epinbu].T
    print 'Refstars filtradas por match: %d' %(rid1.size)
    rid_cat = np.copy(rid1)
    #print 'rid_cat=%d' %(rid_cat.size)
    #Estas son las refstars que aparecen en rojo en el plot: [(rm>mp1)*(rm<mp2)] -> if true & true
    #Falta el constrain relacionado a estar dentro del circulo
    prid,prx1,pry1,prm,prx2,pry2 = np.transpose([rid1,rx1,ry1,rm,rx2,ry2])[(rm>mp1)*(rm<mp2)].T
    #filtro las refstars por magnitud
    rid1,rx1,ry1,rm,rx2,ry2 = np.transpose([rid1,rx1,ry1,rm,rx2,ry2])[(rm>mr1)*(rm<mr2)].T

    print 'Refstars para transformaciones (%.1f<Ks<%.1f): %d' %(mr1,mr2,rid1.size)
    print 'Refstars a incluir en el plot (%.1f<Ks<%.1f): %d' %(mp1,mp2,prid.size)

    #Filtro catalogo por magnitudes
    mask = (m<ma2)*(m>ma1)
    id1,x1,y1,m,x2,y2 = np.transpose([id1,x1,y1,m,x2,y2])[mask].T
    print "Estrellas a analizar (%.1f<Ks<%.1f): %d" %(ma1,ma2,id1.size)
    #como imprimo la ubicacion?
    #print 'id1', id1
    #print np.argwhere (id1==55)
    #print id1==55
    #sys.exit()

    #defino epxy para la funcion que busca refstars cercanas (nbrs)
    epxy = np.transpose([rx2,ry2])

    if local:
        #?esta funcion usa el parametro nrefstars?
        nbrs = NN(n_neighbors=nrefstars, algorithm='auto').fit(epxy)

        #coords transformadas. Se crea arreglo de ceros
        ctx = np.zeros(x1.size)
        cty = np.zeros(y1.size)

        if rad_ext!=0:

            print '\nQuieres %d refstars para transformar cada estrella' %(nrefstars-1)

            dist, nei = nbrs.radius_neighbors(np.transpose([x2,y2]),radius=rad_ext)
            nbors      = np.array([len(d) for d in dist])

            #los nei son tantos arreglos con la posicion (indice) de las refstars
            #como estrellas a analizar
            #si imprimo nei da el indice de las refstars seleccionadas para todas las estrellas
            #print 'nei:', nei
            #print 'nei:', nei.shape
            #print 'nei:', nei[0]
            #print 'refstars a distancia < rad_ext (id) \n', rid1[nei[0]]
            #print 'Numero de refstars encontradas por estrella (r<rad_ext):', nbors

            #los dist son tantos arreglos con las distancias de las refstars como
            #estrellas a analizar
            #?como imprimo las distancias de menor a mayor
            #print 'Distancias de las refstars por estrella:\n', dist[0]
            #sys.exit()

            mednbors, minnbors,maxnbors = np.median(nbors),nbors.min(),nbors.max()
            print'valor min,medio,max de refstars:',minnbors,mednbors,maxnbors
            #sys.exit()

            print '\nElegimos para cada estrella las %d refstars...' %(nrefstars-1)
            if sort_mag:
                print 'mas brillantes'
            else:
                print 'mas cercanas'

            for j in range(len(dist)):
                #Elimina la misma estrella y las interiores a rad_int
                msk     = (dist[j]!=0)*(dist[j]>rad_int)
                dist[j] = dist[j][msk]
                nei[j]  = nei[j][msk]
                #print'j =',j
                #print'nei[j] =', nei
                #print'id', rid1[nei[j]]
                #print'numero de refstars', len(dist[j])
                #print'distancia de refstars\n', (dist[j])
                #Toma las mas brillantes si tengo mas que nrefstars
                if len(dist[j])>nrefstars:
                    if sort_mag:
                        lrm    = rm[nei[j]]
                        #?porque los :?
                        mnei    = np.argsort(lrm)[:nrefstars]
                        dist[j] = dist[j][mnei]
                        nei[j]  = nei[j][mnei]
                #print'numero de refstars', len(dist[j])
                #Toma las mas cercanas si tengo mas que nrefstars
                    else:
                        mnei     = np.argsort(dist[j])[:nrefstars]
                        dist[j] = dist[j][mnei]
                        nei[j]  = nei[j][mnei]
                        #print'numero de refstars', len(dist[j])
                        #con este lo paro cuando termina de analizar j=0
                        #sys.exit()
                        #con este lo paro cuando termina el for
                        #sys.exit()

        #Refstars mas cercanos sin restriccion de distancias (rad_ext=0)
        else:
            dist,nei = nbrs.kneighbors(np.transpose([x2,y2]))
            print '\nBuscando refstars sin restriccion de distancias (rad_ext=0)'
            print 'Estrellas a transformar, refstars por estrella:', nei.shape
            ##recordar: los nei son arreglos con el indice de las refstars
            #print 'refstars mas cercanas', nei
            #print 'refstars mas cercanas (id) \n', id1[nei]
            #print '\ndistancias refstars mas cercanos \n', dist

            #le saca el primero porque puede contener la misma estrella
            nei      = nei[:,1:]
            dist      = dist[:,1:]

        #print '\nrefstars mas cercanas', nei.shape
        #print 'refstars mas cercanas', nei
        #sys.exit()

        #?porque no np.mean(dist)?
        means = []
        for d in dist:
            if len(d)!=0:
                means.append(np.nanmean(d))
        means = np.array(means)
        #means = np.array([np.nanmean(d) for d in dist])
        nbors       = np.array([len(d) for d in dist])
        mednbors, minnbors = np.nanmedian(nbors), np.nanmin(nbors)

        #print '\nctx', ctx.shape
        #print 'ctx', ctx

        #este progressBar muestra las flecha moviendose a medida que hace la regresion para
        #cada estrella del catalogo que se esta analizando
        #x1.size son la cantidad de estrellas del catalogo que se esta analizando (dp de filtar)
        with ProgressBar(x1.size) as bar:
            #para cada estrella del catalogo:
            for k in range(x1.size):
                if len(nei[k]) < 4:
                    ctx[k] = np.nan
                    cty[k] = np.nan
                    #ctx[k] = x1[k] - 888.8
                    #cty[k] = y1[k] - 888.8
                    continue
                    #print '\nERROR: Encontro muy pocas refstars!'
                    #sys.exit()
                #recordar: los nei son arreglos con el indice de las refstars seleccionadas
                #print'\n\nnei[k]', nei[k].shape
                #print'', nei[k]
                #imprime id de refstars (vecinas) a usar
                #print'\nid correspondiente \n', id1[nei[k]]

                #coordenadas originales de las refstars seleccionadas
                coords = np.transpose([rx2,ry2])[nei[k]].T
                #print'\ncoords a transformar', coords.shape
                #print'', coords
                #sys.exit()

                #coordenadas en el sist de referencia de las refstars seleccionadas
                ep1_x = rx1[nei[k]]
                ep1_y = ry1[nei[k]]
                #print'\ncoords de las refstars en reference catalog\n', ep1_x
                #print'', ep1_y
                #sys.exit()

                #queremos resolver rx1=a*rx2+b*ry2+c
                #coords=rx2 ry2, ep1_x/ep1_y son las coord en ref catalog (rx1,ry1)
                poptx, pcovx = curve_fit(linear,coords,ep1_x)
                popty, pcovy = curve_fit(linear,coords,ep1_y)
                #print'\npoptx', poptx
                #print'popty', popty


                        #calcula las coords transformadas (ct):se le dan x2,y2 y usa poptx,popty
                ctx[k] += linear([x2[k],y2[k]],*poptx)
                cty[k] += linear([x2[k],y2[k]],*popty)
                                #print'\nx_refcat:', x1[k]
                                #print'y_refcat:', y1[k]
                #print'\nx_transformed:', ctx[k]
                #print'y_transformed:', cty[k]

                #sys.exit()
                bar.update()

    #Global:No usa las refstars mas cercanas
    #usa todas las x2,y2 para el fit, o sea todas las estrellas del catalogo dp de filtrar
    else:
        print'\nTransformacion Global:usando las %d estrellas filtradas' %(x2.size)
        #?Porque non-linear? ####Use non-linear least squares to fit a function, f, to data.
        poptx, pcovx = curve_fit(linear,[x2,y2],x1)
        popty, pcovy = curve_fit(linear,[x2,y2],y1)
        #print'\nx2', x2
        #print'y2', y2
        #print'\nx1', x1
        #print'y1', y1

        #obtengo los valores del fit ax+by+c
        #print'\npoptx', poptx
        #print'popty', popty
        #matrices de cov
        #print'\npcovx', pcovx
        #print'\npcovy', pcovy

        #calcula las coords transformadas (ct):se le dan x2,y2 y usa poptx,popty
        ctx = linear([x2,y2],*poptx)
        cty = linear([x2,y2],*popty)
        #print'\nx_transformed:\n', ctx
        #print'y_transformed\n', cty

        #cuando se usan? En el plot de output, se dan estos valores para la dist_med,
        #dist_max, media_vec, min_vec
        means = -1
        nbors = -1
        mednbors, minnbors = -1,-1
        #aqui termina Global

    #VALORES PARA TODAS LAS ESTRELLAS
    dx = x1-ctx
    dy = y1-cty

    print '\n        dx, dy'
    print 'num:     %d, %d' % (dx.size,dy.size)
    print 'max:     %8.8f, %8.8f' % (np.nanmax(dx),np.nanmax(dx))
    print 'min:     %8.8f, %8.8f' % (np.nanmin(dy),np.nanmin(dy))
    print 'mean:    %8.8f, %8.8f' % (np.nanmean(dx),np.nanmean(dy))
    print 'median:  %8.8f, %8.8f' % (np.nanmedian(dx),np.nanmedian(dy))
    print 'std:     %8.8f, %8.8f' % (np.nanstd(dx),np.nanstd(dy))

    #Terminadas las transformaciones (Local o Global)
    #se guarda id, delta_x y delta_y
    ctx[np.isnan(ctx)] = x1[np.isnan(ctx)] - 888.8
    cty[np.isnan(cty)] = y1[np.isnan(cty)] - 888.8
    data = np.transpose([id1,x1-ctx,y1-cty])
    #print'\ndata', data.shape
    #print'\nid delta_x delta_y \n', data
    #sys.exit()
    makedir('PMs')
    np.savetxt('./PMs/PM_%03d.dat' % nro_epoca[i], data, header='ID DX DY', fmt='%d %f %f')

    #PLOT OUTPUT.PSF

    #print 'Estrellas analizadas (id1)', id1.size
    #print 'Refstars usadas (rid1)', rid1.size
    #print 'Refstars en mi catalogo (rid_cat)', rid_cat.size

    #ref son solo true y false, tantos como elementos en id1(catalogo). El 1ero manda
    #ref  = np.in1d(id1,rid_cat) ###este es original de Felipe
    #m_ref son solo true y false, tantos como rm o m, segun sea el caso. Los true son los que
    #estan en el rango de mp1 y mp2
    #m_ref  = (m>mp1)*(m<mp2) ###este es original de Felipe
    #r_ref  = (np.sqrt((x1-x0)**2 + (y1-y0)**2) < rad_ref) ###este es original de Felipe
    #print'\nref:', ref.shape
    #print'm_ref:', m_ref.shape
    #print'r_ref:', r_ref.shape

    #rplot   = np.size((x1)[ref])
    #print '\nrpoints_1:', rplot
    #rplot   = np.size((x1)[ref*m_ref*r_ref])
    #print 'rpoints_2:', rplot
    #bplot   = np.size((x1)[~ref])
    #print 'bpoints:', bplot

    ref  = np.in1d(id1,rid1)
    m_ref  = (m>mp1)*(m<mp2)
    r_ref  = (np.sqrt((x1-x0)**2 + (y1-y0)**2) < rad_ref) ###este es original de Felipe
    #print'\nref:', ref.shape
    #print'm_ref:', m_ref.shape
    #rplot   = np.size((x1)[ref])
    #print '\nrpoints_1:', rplot
    #rplot   = np.size((x1)[ref*m_ref*r_ref])
    #print 'rpoints_2:', rplot
    #bplot   = np.size((x1)[~ref])
    #print 'bpoints:', bplot
    #sys.exit()

    #####PTOS AZULES (CATALOGO-REFSTARS)
    a.scatter((x1-ctx)[~ref],(y1-cty)[~ref],s=1,rasterized=True,edgecolor='',color='#0055FF',lw=.5)
    #este lo puse yo para plotear circulos
    #a.scatter((x1-ctx)[~ref],(y1-cty)[~ref],s=1, facecolors='none', edgecolors='b',lw=.5)

    #####PTOS ROJOS (REFSTARS)
    a.scatter((x1-ctx)[ref*m_ref*r_ref],(y1-cty)[ref*m_ref*r_ref],s=1.25,rasterized=True,edgecolor='',color='#FF5500',lw=.5)

    #TEXTO INFERIOR EN PLOT output.pdf
    a.text(.05,.05,u'$S = %f$\n$E = %s$\n$Nr/Nb = %d/%d$\n$Med_d, Max_d = %.3f, %.3f$\n$Med_n,Min_n = %d,%d$'%(se[i],el[i],(ref*m_ref*r_ref).sum(),(~ref).sum(),np.nanmedian(means),np.nanmax(means),mednbors, minnbors),transform = a.transAxes,alpha=.66,fontsize=10)
    a.text(.85,.9,u'$%d$' % nro_epoca[i],transform = a.transAxes,alpha=.66,fontsize=14)

    #sys.exit()

    #Estadistica estrellas azules a incluir en plot output.psf
    xmean_b = np.nanmean((x1-ctx)[~ref])
    ymean_b = np.nanmean((y1-cty)[~ref])
    xstd_b    = np.nanstd((x1-ctx)[~ref])
    ystd_b    = np.nanstd((y1-cty)[~ref])

    #Estadistica estrellas rojas (refstars) a incluir en plot output.psf
    #se cae si uso los .gc y rad_ref o x0,y0 referido a los archivos.mat
    xmean_r = np.nanmean((x1-ctx)[ref*m_ref*r_ref])
    ymean_r = np.nanmean((y1-cty)[ref*m_ref*r_ref])
    xstd_r    = np.nanstd((x1-ctx)[ref*m_ref*r_ref])
    ystd_r    = np.nanstd((y1-cty)[ref*m_ref*r_ref])
    #print 'xmean_r', xmean_r

    #?Porque definir esto?
    meansx_fie[i] += xmean_b
    meansy_fie[i] += ymean_b
    stdx_fie[i]   += xstd_b
    stdy_fie[i]   += ystd_b

    meansx_ref[i] += xmean_r
    meansy_ref[i] += ymean_r
    stdx_ref[i]   += xstd_r
    stdy_ref[i]   += ystd_r

    #?como se imprime la epoca?
    #TEXTO SUPERIOR EN PLOT output.pdf
    a.text(.05,.85,u'M = $%.4f / %.4f$\nS = $%.3f / %.3f$' % (xmean_r,ymean_r,xstd_r,ystd_r),color='#FF5500',transform=a.transAxes,fontsize=10)
    a.text(.05,.73,u'M = $%.4f / %.4f$\nS = $%.3f / %.3f$' % (xmean_b,ymean_b,xstd_b,ystd_b),color='#0055FF',transform=a.transAxes,fontsize=10)
    a.set_xlim(-lim,lim)
    a.set_ylim(-lim,lim)
    a.set_aspect('equal')

    #Cambia tamano de la fuente de los ticks
    #?en cual plot se usa esta instruccion?
    a.tick_params(axis='both', which='major', labelsize=8)

    #PLOT OUTPUT_DEL_XY.PSF
    ad[2*i].scatter(x1[~ref],dx[~ref],s=1,rasterized=True,lw=0,color='#0055FF')
    ad[2*i].scatter(x1[ref*m_ref],dx[ref*m_ref],s=1,rasterized=True,lw=0,color='#FF5500')
    ##esto es nuevo
    ad[2*i].tick_params(axis='both', which='major', labelsize=8)
    ad[2*i].set_xlim(x1.min()-100,x1.max()+100)
    ad[2*i].set_ylim(-1,1)

    ad[2*i+1].scatter(y1[~ref],dy[~ref],s=1,rasterized=True,lw=0,color='#0055FF')
    ad[2*i+1].scatter(y1[ref*m_ref],dy[ref*m_ref],s=1,rasterized=True,lw=0,color='#FF5500')
    ##esto es nuevo
    ad[2*i+1].tick_params(axis='both', which='major', labelsize=8)
    ad[2*i+1].set_xlim(y1.min()-100, y1.max()+100)
    ad[2*i].set_ylim(-1,1)
    ad[2*i+1].text(.9,.8,'%d' % nro_epoca[i], transform=ad[2*i+1].transAxes, fontsize=8)

    if i==2:
        print '\nGuardando plot de primeras 3 epocas...\n'
        try:
            muygrande = False
            fig.savefig(output+'.pdf',dpi=200)
            if plot_del_xy:
                #?esto es nuevo...a que sirve?
                fig_delta.tight_layout()
                fig_delta.savefig(output+'_del_xy.pdf',dpi=200)
        except:
            print '\nAdvertencia: El plot final es demasiado grande, por lo que no se guardara!'
            muygrande = True

#Aca termina el for para todos los catalogos de input

print '\nGuardando plot final...'

#HEADER PLOT OUTPUT.PSF
fig.suptitle('Refstars:%3d; Radios:%3d,%4d; Mag stars:%2.1f,%2.1f; Mag refstars:%2.1f,%2.1f; Mag ref_plot:%2.1f,%2.1f' % (nrefstars, rad_int, rad_ext, ma1, ma2, mr1, mr2, mp1, mp2))
fig.tight_layout()
fig.subplots_adjust(top=0.95)
if not muygrande:
    fig.savefig(output+'.pdf',dpi=200)
    fig_delta.savefig(output+'_del_xy.pdf',dpi=200)
else:
    plt.show()

yrs = (yr-yr[0])/365.25
eff_yrs = yrs[nro_epoca-1]
#print 'yrs:', yrs
#print 'eff_yrs:', eff_yrs
#sys.exit()

#PLOT DELTAS VS TIEMPO
if plot_del_ep:
    if len(nro_epoca) < 2:
        print "Hay solo un catalogo. El plot de delta vs tiempo no se genera."
        sys.exit()

    fig, ax = plt.subplots(nrows=2,ncols=2,figsize=[6*2,2*2])

    ax[0,0].errorbar(eff_yrs, meansx_ref, yerr=stdx_ref ,fmt='.', c='#FF5500', ms=13, rasterized=True)
    ax[0,1].errorbar(eff_yrs, meansy_ref, yerr=stdy_ref, fmt='.', c='#FF5500', ms=13, rasterized=True)

    ax[1,0].errorbar(eff_yrs, meansx_fie, yerr=stdx_fie, fmt='.', c='#0055FF', ms=13, rasterized=True)
    ax[1,1].errorbar(eff_yrs, meansy_fie, yerr=stdy_fie, fmt='.', c='#0055FF', ms=13, rasterized=True)

    poptrxc, pcovrxc = curve_fit(recta, eff_yrs, meansx_ref, sigma=stdx_ref, absolute_sigma=True)
    poptryc, pcovryc = curve_fit(recta, eff_yrs, meansy_ref, sigma=stdy_ref, absolute_sigma=True)

    poptrxf, pcovrxf = curve_fit(recta, eff_yrs, meansx_fie, sigma=stdx_fie, absolute_sigma=True)
    poptryf, pcovryf = curve_fit(recta, eff_yrs, meansy_fie, sigma=stdy_fie, absolute_sigma=True)

    ax[0,0].plot(eff_yrs, recta(eff_yrs, *poptrxc))
    ax[0,0].text(.6,.8,r'$y = %fx + %f$' % tuple(poptrxc),transform=ax[0,0].transAxes)
    ax[0,1].plot(eff_yrs, recta(eff_yrs, *poptryc))
    ax[0,1].text(.6,.8,r'$y = %fx + %f$' % tuple(poptryc),transform=ax[0,1].transAxes)

    ax[1,0].plot(eff_yrs, recta(eff_yrs, *poptrxf))
    ax[1,0].text(.6,.8,r'$y = %fx + %f$' % tuple(poptrxf),transform=ax[1,0].transAxes)
    ax[1,1].plot(eff_yrs, recta(eff_yrs, *poptryf))
    ax[1,1].text(.6,.8,r'$y = %fx + %f$' % tuple(poptryf),transform=ax[1,1].transAxes)

    diferencias = np.zeros(4)
    for i,a in enumerate(np.ravel(ax)):
        rmin, rmax = a.get_ylim()
        diferencias[i] += np.abs(rmax-rmin)

    max_dif = np.nanmax(diferencias)
    add_dif = (max_dif - diferencias) / 2.

    for i,a in enumerate(np.ravel(ax)):
        rmin, rmax = a.get_ylim()
        a.set_ylim(rmin - add_dif[i], rmax + add_dif[i])

    fig.suptitle('\nRefstars:%3d; Radios:%3d,%4d; Mags:%2.1f,%2.1f; Mag refstars:%2.1f,%2.1f' % (nrefstars, rad_int, rad_ext, ma1, ma2, mr1, mr2) )
    fig.tight_layout()
    fig.subplots_adjust(top=0.85)
    fig.savefig(output+'_del_ep.pdf',dpi=200)
