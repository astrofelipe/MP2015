import sys
import argparse
import os
import subprocess
import numpy as np
import pm_funcs
from astropy.utils.console import color_print
from scipy.optimize import curve_fit
from sklearn.neighbors import KernelDensity

#PARAMETROS
radio, itera, output, refer, nframes, min_ep, max_err = pm_funcs.get_script()

print 'Iniciando script con...'
print 'linea de comando: %s' % ' '.join(sys.argv)
print 'radio: %d' % radio
print 'itera: %d' % itera
print

#ARGUMENTOS
parser = argparse.ArgumentParser(description='Script PM VVV')
parser.add_argument('<Input List>', help='Lista con inputs (para tlineal)')
parser.add_argument('<Ref Catalog>', help='Catalogo de referencia (usado por PM_1a1)')
parser.add_argument('-c', '--continua', type=int, default=None, help='Realiza n iteraciones partiendo desde la ultima hecha anteriormente')
parser.add_argument('-p', '--peak', type=int, help='Usa el peak de los PM para calcular nuevas refstars')
parser.add_argument('-r', '--refstars', action='store_true', help='Usa VPD + refstars0, sino solo VPD')

args = parser.parse_args()

#FUNCIONES Y OTROS PARAMETROS
inputs        = vars(args)['<Input List>']#sys.argv[1]
ref_cat       = vars(args)['<Ref Catalog>']#sys.argv[2]
stilts_folder = os.path.dirname(os.path.realpath(__file__))

#continua = False
#if (len(sys.argv) > 3) and (sys.argv[3]=='-C'):
#    continua = True
#continua = args.continua
if args.continua is not None:
    continua = True
    itera    = args.continua
else:
    continua = False

def makedir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def gaussian(x, amp, mu, sig):
    if (amp < 0) or (sig < 0):
        return np.inf
    return amp * np.exp(-(x-mu)**2 / (2*sig**2))

def two_gaussian(x, amp1, mu1, sig1,
                    amp2, mu2, sig2):
    return (gaussian(x, amp1, mu1, sig1) +
            gaussian(x, amp2, mu2, sig2))

def three_gaussian(x, amp1, mu1, sig1,
                      amp2, mu2, sig2,
                      amp3, mu3, sig3):
    return (gaussian(x, amp1, mu1, sig1) +
    gaussian(x, amp2, mu2, sig2) +
    gaussian(x, amp3, mu3, sig3))

if args.peak == 1:
    gf = gaussian
if args.peak == 2:
    gf = two_gaussian
if args.peak == 3:
    gf = three_gaussian



#PIPELINE
#output = subprocess.check_output('grep "PDF de Output" %s/tlineal_1a1.py' % stilts_folder, shell=True)
#output = output.split("'")[1]

#refstars = subprocess.check_output('grep "Catalogo con las estrellas de referencia" %s/tlineal_1a1.py' % stilts_folder, shell=True)
#refstars = refstars.split("'")[1]

#nframes = subprocess.check_output('grep "Numero minimo de epocas en que debe estar la estrella" %s/pm_1a1.py' % stilts_folder, shell=True)
#nframes = int(nframes.split(' ')[2])

#min_ep = subprocess.check_output('grep "min_ep =" %s/VPDHmag.py' % stilts_folder, shell=True)
#min_ep = int(min_ep.split(' ')[-1])
refstars = refer
print '\nOutput en tlineal_1a1.py: %s' % output
print 'Archivo de refstars: %s' % refer
print 'nframes en pm_1a1.py: %d' % nframes
print 'min_ep en VPDHmag.py: %d' % min_ep
if continua:
    print 'Continuando desde la ultima iteracion'

nro_files = np.genfromtxt(inputs, unpack=True, dtype='string')
if not nro_files.shape:
    nro_files = np.atleast_1d(nro_files)
nro_files = nro_files.size

if nro_files < nframes:
    print '\nNumero de archivos en el input es menor que nframes!'
    sys.exit()

if nro_files < min_ep:
    print '\nNumero de archivos en el input es menor que min_ep!'
    sys.exit()

if not os.path.isfile('refstars0.gc'):
    print '\nArchivo refstars0.gc no encontrado!'
    sys.exit()

idsr = np.genfromtxt('refstars0.gc', unpack=True, usecols=(0,))

if not continua:
    if os.path.isfile('zelimchisha'):
        print '\nEliminando estrellas de zelimchisha'
        rej_ids = np.genfromtxt('zelimchisha', unpack=True, usecols=(0,))
        id_mask = ~np.in1d(idsr, rej_ids)
        print '\tEncontradas %d estrellas' % (~id_mask).sum()
    else:
        id_mask = np.ones(len(idsr)).astype(bool)

    data = np.genfromtxt('refstars0.gc')
    data = data[id_mask]

    print 'Guardando refstars para primera iteracion'
    np.savetxt(refstars, data, fmt='%f')

    last_idx = 0

    #for i in range(1,itera+1):
    #    if os.path.exists('iter_%d' % i):
    subprocess.call('rm -r iter_*', shell=True)

if continua:
    import glob
    iter_dirs = np.sort(glob.glob('iter_*'))

    if len(iter_dirs) == 0:
        print 'No hay carpetas iter!'
        sys.exit(1)
    last_idx  = int(iter_dirs[-1].split('_')[-1])

    if not os.path.isfile('iter_%d/PM_final.dat' % last_idx):
        print '\nNo se encontro PM_final.dat en iter_%d!' % last_idx
        sys.exit(1)

    print '\nGenerando nuevo archivo de refstars a partir de la ultima iteracion'
    ids, pmx, pmy, nf, pmxe, pmye = np.genfromtxt('iter_%d/PM_final.dat' % (last_idx), unpack=True, usecols=(0,3,4,6,8,9))
    #Filtro por errores
    pme = np.sqrt(pmxe**2 + pmye**2)
    pm1 = (pmx**2 + pmy**2)**0.5 < 30

    #Filtro por zelimchisha
    if os.path.isfile('zelimchisha'):
        rej_ids = np.genfromtxt('zelimchisha', unpack=True, usecols=(0,))
        id_mask = ~np.in1d(ids, rej_ids)
    else:
        id_mask = np.ones(len(ids)).astype(bool)
    #Centrar PMs
    if args.peak != None:
        print '\tCalculando peak para centrar refstars'
        from sklearn import mixture
        XG    = np.transpose([pmx, pmy])[(nf >= nframes) * (pme <= max_err) * (id_mask) * (pm1)]
        gmm   = mixture.GMM(n_components=args.peak, covariance_type='full').fit(XG)
        xg,yg = np.transpose(gmm.means_)

        #Centros x
        kde = KernelDensity(kernel='gaussian').fit(XG[:,0][:, np.newaxis])
        xx  = np.arange(-15, 15, 0.05)
        yy  = np.exp(kde.score_samples(xx[:,np.newaxis]))
        a0  = np.exp(kde.score_samples(xg[:,np.newaxis]))
        if args.peak == 3:
            p0  = [a0[0]/2.0, xg[0], 3, a0[1]/2.0, xg[1], 3, a0[2]/2.0, xg[2], 3]
        elif args.peak == 2:
            p0  = [a0[0]/2.0, xg[0], 3, a0[1]/2.0, xg[1], 3]
        elif args.peak == 1:
            p0  = [a0[0]/2.0, xg[0], 3]

        popt, pcov = curve_fit(gf, xx, yy, p0=p0, maxfev=10000)
        x0g = popt[1::3]
        x0e = np.sqrt(np.diag(pcov)[1::3])

        #Centros y
        kde = KernelDensity(kernel='gaussian').fit(XG[:,1][:, np.newaxis])
        yy  = np.exp(kde.score_samples(xx[:,np.newaxis]))
        a0  = np.exp(kde.score_samples(yg[:,np.newaxis]))
        if args.peak == 3:
            p0  = [a0[0]/2.0, yg[0], 3, a0[1]/2.0, yg[1], 3, a0[2]/2.0, yg[2], 3]
        elif args.peak == 2:
            p0  = [a0[0]/2.0, yg[0], 3, a0[1]/2.0, yg[1], 3]
        elif args.peak == 1:
            p0  = [a0[0]/2.0, yg[0], 3]

        popt, pcov = curve_fit(gf, xx, yy, p0=p0, maxfev=10000)
        y0g = popt[1::3]
        y0e = np.sqrt(np.diag(pcov)[1::3])

        cen = np.transpose([x0g, y0g])
        idx = np.argmin(np.sum(cen**2, axis=1))

        if np.sum(cen[idx]**2) > 3:
            x0,y0 = gmm.means_[idx]
            print '\tCentros muy grandes, volviendo a los guess'
        else:
            x0,y0 = cen[idx]
        print '\t', x0, y0
    else:
        x0,y0 = 0,0
    #Filtro por radio en PM
    pmr = np.sqrt((pmx-x0)**2 + (pmy-y0)**2)

    if args.refstars:
        pm0 = np.in1d(ids, idsr)
    else:
        pm0 = True

    mask = (pmr <= radio) * (nf >= nframes) * (pme <= max_err) * (id_mask) * (pm0)

    data = np.genfromtxt('iter_%d/PM_final.dat' % (last_idx))
    data = data[mask]

    fmt = '%d %.6f %.6f %.6f %.6f %.3f %d %d %.6f %.6f %.0f %.2f'
    hdr = 'ID RA DEC PM_X PM_Y MAG_K NFRAMES CFRAMES PMXE PMYE NEI NEI_STD'

    np.savetxt(refstars, data, fmt=fmt, header=hdr)

for i in range(itera):
    color_print('\nComenzando iteracion: %d' % (last_idx+i+1), 'lightcyan')

    #Crea carpeta para guardar los outputs
    makedir('iter_%d' % (last_idx+i+1))

    color_print('\tEjecutando tlineal_1a1.py', 'cyan')
    subprocess.call('python -u %s/tlineal_1a1.py %s' % (stilts_folder, inputs), shell=True)

    color_print('\tEjecutando pm_1a1.py', 'cyan')
    subprocess.call('python %s/pm_1a1.py %s' % (stilts_folder, ref_cat), shell=True)

    color_print('\tEjecutando VPDHmag.py', 'cyan')
    subprocess.call('python %s/VPDHmag.py' % stilts_folder, shell=True)

    color_print('\tMoviendo archivos', 'cyan')
    subprocess.call('mv %s*.png iter_%d' % (output, (last_idx+i+1)), shell=True)
    subprocess.call('mv %s iter_%d' % (refstars, (last_idx+i+1)), shell=True)
    makedir('iter_%d/PMs' % (last_idx+i+1))
    subprocess.call('mv PM*.h5 iter_%d' % (last_idx+i+1), shell=True)
    subprocess.call('mv PM_final.dat iter_%d' % (last_idx+i+1), shell=True)
    subprocess.call('mv PMs/* iter_%d/PMs' % (last_idx+i+1), shell=True)
    subprocess.call('mv VPD*.png iter_%d' % (last_idx+i+1), shell=True)

    color_print('\tGenerando nuevo archivo de refstars', 'cyan')
    ids, pmx, pmy, nf, pmex, pmey = np.genfromtxt('iter_%d/PM_final.dat' % (last_idx+i+1), unpack=True, usecols=(0,3,4,6,8,9))
    pme = np.sqrt(pmex**2 + pmey**2)
    pm1 = (pmx**2 + pmy**2)**0.5 < 30


    if os.path.isfile('zelimchisha'):
        rej_ids = np.genfromtxt('zelimchisha', unpack=True, usecols=(0,))
        id_mask = ~np.in1d(ids, rej_ids)
    else:
        id_mask = np.ones(len(ids)).astype(bool)
    #Centrar PMs
    if args.peak != None:
        print '\tCalculando peak para centrar refstars'
        from sklearn import mixture
        XG    = np.transpose([pmx, pmy])[(nf >= nframes) * (pme <= max_err) * (id_mask) * (pm1)]
        gmm   = mixture.GMM(n_components=args.peak, covariance_type='full').fit(XG)
        xg,yg = np.transpose(gmm.means_)

        #Centros x
        kde = KernelDensity(kernel='gaussian').fit(XG[:,0][:, np.newaxis])
        xx  = np.arange(-15, 15, 0.05)
        yy  = np.exp(kde.score_samples(xx[:,np.newaxis]))
        a0  = np.exp(kde.score_samples(xg[:,np.newaxis]))
        if args.peak == 3:
            p0  = [a0[0]/2.0, xg[0], 3, a0[1]/2.0, xg[1], 3, a0[2]/2.0, xg[2], 3]
        elif args.peak == 2:
            p0  = [a0[0]/2.0, xg[0], 3, a0[1]/2.0, xg[1], 3]
        elif args.peak == 1:
            p0  = [a0[0]/2.0, xg[0], 3]

        popt, pcov = curve_fit(gf, xx, yy, p0=p0, maxfev=10000)
        x0g = popt[1::3]
        x0e = np.sqrt(np.diag(pcov)[1::3])

        #Centros y
        kde = KernelDensity(kernel='gaussian').fit(XG[:,1][:, np.newaxis])
        yy  = np.exp(kde.score_samples(xx[:,np.newaxis]))
        a0  = np.exp(kde.score_samples(yg[:,np.newaxis]))
        if args.peak == 3:
            p0  = [a0[0]/2.0, yg[0], 3, a0[1]/2.0, yg[1], 3, a0[2]/2.0, yg[2], 3]
        elif args.peak == 2:
            p0  = [a0[0]/2.0, yg[0], 3, a0[1]/2.0, yg[1], 3]
        elif args.peak == 1:
            p0  = [a0[0]/2.0, yg[0], 3]

        popt, pcov = curve_fit(gf, xx, yy, p0=p0, maxfev=10000)
        y0g = popt[1::3]
        y0e = np.sqrt(np.diag(pcov)[1::3])

        cen = np.transpose([x0g, y0g])
        idx = np.argmin(np.sum(cen**2, axis=1))

        if np.sum(cen[idx]**2) > 3:
            x0,y0 = gmm.means_[idx]
            print '\tCentros muy grandes, volviendo a los guess'
        else:
            x0,y0 = cen[idx]
        print '\t', x0, y0
    else:
        x0,y0 = 0,0

    pmr = np.sqrt((pmx-x0)**2 + (pmy-y0)**2)

    if args.refstars:
        pm0 = np.in1d(ids, idsr)
    else:
        pm0 = True

    mask = (pmr <= radio) * (nf >= nframes) * (pme <= max_err) * (id_mask) * (pm0)

    data = np.genfromtxt('iter_%d/PM_final.dat' % (last_idx+i+1))
    data = data[mask]

    fmt = '%d %.6f %.6f %.6f %.6f %.3f %d %d %.6f %.6f %.0f %.2f'
    hdr = 'ID RA DEC PM_X PM_Y MAG_K NFRAMES CFRAMES PMXE PMYE NEI NEI_STD'

    np.savetxt(refstars, data, fmt=fmt, header=hdr)

subprocess.call('rm -r PMs', shell=True)
subprocess.call('rm %s' % refstars, shell=True)

print '\nDone!'
