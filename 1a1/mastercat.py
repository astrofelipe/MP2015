import numpy as np
from scipy.optimize import curve_fit
import sys
from joblib import Parallel, delayed
from sklearn.neighbors import NearestNeighbors as NN

min_epochs = 5      #Numero minimo de epocas para considerar la estrella
min_mag    = 11     #Magnitud minima para estrellas a transformar
max_mag    = 14     #Magnitud maxima...
max_err    = .05    #Error maximo a considerar
iteraciones = 5
iteracion2 = 'global'

nrefstars = 51 #Numero de vecinos locales (contando a la estrella propia) si iteracion2==global

masterst   = sys.argv[1]

#FUNCIONES
def linear(coords, a, b, c):
    x, y = coords
    return a*x + b*y + c

#PIPELINE
data = np.genfromtxt(masterst, unpack=True, delimiter=',')
data = data.T

ids = data[:, 0::7]

nro_ep = ids.shape[1]
found  = np.sum(np.isfinite(ids), axis=1)
rej    = found >= min_epochs

print 'Numero de epocas: %d' % nro_ep
print 'Numero total de estrellas: %d' % len(ids)
print 'Numero de estrellas en %d epocas: %d' % (min_epochs, np.sum(rej))

ids = data[:, 0::7][rej]
ras = data[:, 1::7][rej]
des = data[:, 2::7][rej]
xs  = data[:, 3::7][rej]
ys  = data[:, 4::7][rej]
ms  = data[:, 5::7][rej]
es  = data[:, 6::7][rej]

def transformacion(ep):
    magcon = (ms[:, 0] > min_mag) * (ms[:, 0] < max_mag) * (es[:, 0] < max_err)
    common = np.isfinite(ids[:, 0]) * np.isfinite(ids[:, ep])
    print 'Numero de estrellas utilizadas en %d: %d' % (ep, np.sum(magcon*common))

    x1 = xs[:, 0]
    x2 = xs[:, ep]
    y1 = ys[:, 0]
    y2 = ys[:, ep]

    xx1 = x1[common*magcon]
    xx2 = x2[common*magcon]
    yy1 = y1[common*magcon]
    yy2 = y2[common*magcon]

    m1 = ms[:, 0][common*magcon]
    m2 = ms[:, ep][common*magcon]

    mdm = np.mean(m1-m2)
    mm2 = ms[:, ep] + mdm

    poptx, pcovx = curve_fit(linear, [xx2, yy2], xx1)
    popty, pcovy = curve_fit(linear, [xx2, yy2], yy1)

    tx = linear([x2,y2], *poptx)
    ty = linear([x2,y2], *popty)

    return tx, ty, mm2

xx = np.empty_like(xs)
yy = np.empty_like(ys)
mm = np.empty_like(ms)
ee = np.copy(es)

xx[:, 0] = xs[:, 0]
yy[:, 0] = ys[:, 0]
mm[:, 0] = ms[:, 0]

print '\nIteracion 1'
for j in range(1, nro_ep):
    tx, ty, mm2 = transformacion(j)

    xx[:,j] = tx
    yy[:,j] = ty
    mm[:,j] = mm2

xs[:,0] = np.nanmean(xx, axis=1)
ys[:,0] = np.nanmean(yy, axis=1)

msmask  = np.ma.array(mm, mask=np.isnan(mm))
ms[:,0] = np.ma.average(msmask, axis=1, weights=1.0/ee)

if iteracion2=='global':
    for i in range(iteraciones-1):
        print '\nIteracion: %d' % (i+2)
        for j in range(1, nro_ep):
            tx, ty, mm2 = transformacion(j)

            xx[:,j] = tx
            yy[:,j] = ty
            mm[:,j] = mm2

        print xs[:,0]
        xs[:,0] = np.nanmean(xx, axis=1)
        print xs[:,0]
        ys[:,0] = np.nanmean(yy, axis=1)
        msmask  = np.ma.array(mm, mask=np.isnan(mm))
        ms[:,0] = np.ma.average(msmask, axis=1, weights=1.0/ee)

if iteracion2=='local':
    for i in range(iteraciones-1):
        print '\nIteracion: %d' % (i+2)
        for j in range(1, nro_ep):
            magcon = (ms[:, 0] > min_mag) * (ms[:, 0] < max_mag) * (es[:, 0] < max_err)
            common = np.isfinite(ids[:, 0]) * np.isfinite(ids[:, j])
            in2    = np.isfinite(ids[:, j])

            x1 = xs[:, 0]
            x2 = xs[:, j]
            y1 = ys[:, 0]
            y2 = ys[:, j]

            xx1 = x1[common*magcon]
            xx2 = x2[common*magcon]
            yy1 = y1[common*magcon]
            yy2 = y2[common*magcon]

            m1 = ms[:, 0][common*magcon]
            m2 = ms[:, j][common*magcon]

            mdm = np.mean(m1-m2)
            mm2 = ms[:, j] + mdm

            #Busca vecinos mas cercanos (en 2)
            epxy = np.transpose([xx2, yy2])
            nbrs = NN(n_neighbors=nrefstars, algorithm='auto').fit(epxy)

            dist, nei = nbrs.kneighbors(np.transpose([x2[in2],y2[in2]]))

            def tl(ii):
                nnei = nei[ii][1:]

                xl1, yl1 = xx1[nnei], yy1[nnei]
                xl2, yl2 = xx2[nnei], yy2[nnei]

                poptx, pcovx = curve_fit(linear, [xl2, yl2], xl1)
                popty, pcovy = curve_fit(linear, [xl2, yl2], yl1)

                tx = linear([x2[in2][ii], y2[in2][ii]], *poptx)
                ty = linear([x2[in2][ii], y2[in2][ii]], *popty)

                return tx, ty

            print 'Procesando epoca %d/%d' % (j, nro_ep-1)
            rr = Parallel(n_jobs=4, verbose=0)(delayed(tl)(ii) for ii in range(len(x2[in2])))
            tx, ty = np.array(zip(*rr))

            xx[:,j][in2] = tx
            yy[:,j][in2] = ty
            mm[:,j]      = mm2
        xs[:,0] = np.nanmean(xx, axis=1)
        ys[:,0] = np.nanmean(yy, axis=1)
        msmask  = np.ma.array(mm, mask=np.isnan(mm))
        ms[:,0] = np.ma.average(msmask, axis=1, weights=1.0/ee)

unique, counts = np.unique(ids[np.isfinite(ids)], return_counts=True)

if np.all(counts==1):
    for i in range(1,ids.shape[1]):
        idx = np.isnan(ids[:,0])
        ids[:,0][idx] = ids[:,i][idx]

else:
    idx = np.isnan(ids[:, 0])
    nid = np.arange(np.sum(idx)) + 1e6
    ids[:, 0][idx] = nid

#Asigna RA DEC
idx   = np.isnan(ras[:, 0])
hayra = np.isfinite(ras)
newra = np.array([ras[i][hayra[i]][0] for i in xrange(len(ras))])
newde = np.array([des[i][hayra[i]][0] for i in xrange(len(des))])

ras[:,0] = newra
des[:,0] = newde

#ERROR FINAL
eemask  = np.ma.array(ee, mask=np.isnan(ee))
es[:,0] = np.sqrt(np.sum(np.square(eemask), axis=1))

final_data = np.array([ids, ras, des, xs, ys, ms, es])[:,:,0]
header     = 'ID RA DEC X Y MAG MAG_ERR'
fmt        = '%.0f %.7f %.7f %.3f %.3f %.3f %.3f'

np.savetxt('MASTERCAT.dat', final_data.T, header=header, fmt=fmt)
