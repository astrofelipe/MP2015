from configparser import ConfigParser
from astropy.utils.console import color_print
import os
import numpy as np

params = ConfigParser()

if os.path.isfile('zparams_pm.py'):
    params.read('zparams_pm.py')
else:
    folder = os.path.dirname(os.path.realpath(__file__))
    params.read('%s/zparams_pm.py' % folder)
    print ''
    color_print('###############################################', 'yellow')
    color_print('#', 'yellow', '                                             ', 'lightcyan', '#', 'yellow')
    color_print('#', 'yellow', '  zparams_pm.py no encontrado en la carpeta  ', 'white', '#', 'yellow')
    color_print('#', 'yellow', '\033[1m' + '      Utilizando parametros por defecto      ', 'lightred', '#', 'yellow')
    color_print('#', 'yellow', '                                             ', 'lightcyan', '#', 'yellow')
    color_print('###############################################', 'yellow')
    print ''

#Regresion Lineal
def linear_regression(x, y, w):
    A     = np.vander(x,2)
    W     = np.diag(w)
    ATWA  = np.dot(A.T, np.dot(W, A))
    #ATA   = np.dot(A.T, A / yerr[:, np.newaxis]**2)
    sig_w = np.linalg.inv(ATWA)
    mu_w  = np.linalg.solve(ATWA, np.dot(A.T, np.dot(W, y)))

    return mu_w, sig_w


def get_script():
    radio, itera, max_err = np.array(params['SCRIPT.PY'].values())
    radio = float(radio)
    itera = int(itera)
    max_err = float(max_err)
    output  = str(params['TLINEAL_1a1.PY']['output'])
    refer   = str(params['TLINEAL_1a1.PY']['refer'])
    nframes = int(params['PM_1a1.PY']['nframes'])
    min_ep  = int(params['VPDHmag.PY']['min_ep'])

    return radio, itera, output, refer, nframes, min_ep, max_err

def get_VPDHmag():
    limplotvp, magl, magh, delta, min_ep, min_nei, sigma_err = np.array(params['VPDHmag.PY'].values()).astype('float')
    mags = magl, magh

    return limplotvp, mags, delta, min_ep, min_nei, sigma_err

def get_tlineal():
    nrefstars_tl, min_nei, rad_int, rad_ext, output, refer, sort_mag, \
    local, ma1, ma2, mr1, mr2, mp1, mp2, rad_ref, x0, y0, lim, plot_ep, plot_del_ep, plot_del_xy, nprocs_tl = params['TLINEAL_1a1.PY'].values()

    nrefstars_tl, min_nei, rad_int, rad_ext, ma1, ma2, mr1, mr2, mp1, mp2, rad_ref, x0, y0, lim, nprocs_tl = \
    np.array([nrefstars_tl, min_nei, rad_int, rad_ext, ma1, ma2, mr1, mr2, mp1, mp2, rad_ref, x0, y0, lim, nprocs_tl]).astype(float)

    refer  = str(refer)
    output = str(output)
    nprocs_tl = int(nprocs_tl)

    sort_mag    = str(sort_mag) == 'True'
    local       = str(local) == 'True'
    plot_ep     = str(plot_ep) == 'True'
    plot_del_ep = str(plot_del_ep) == 'True'
    plot_del_xy = str(plot_del_xy) == 'True'

    return nrefstars_tl, min_nei, rad_int, rad_ext, output, refer, sort_mag, \
    local, ma1, ma2, mr1, mr2, mp1, mp2, rad_ref, x0, y0, lim, plot_ep, plot_del_ep, plot_del_xy, nprocs_tl

def get_CMD():

    cmd_modo, match, col1, col2, mag1, mag2, cmd_pdf = params['CMD.PY'].values()
    col1, col2, mag1, mag2 = np.array([col1, col2, mag1, mag2]).astype(float)
    cmd_modo = str(cmd_modo)
    match    = str(match)
    cmd_pdf  = str(cmd_pdf) == 'True'

    return cmd_modo, match, col1, col2, mag1, mag2, cmd_pdf

def get_master_stilts():
    match_tol = float(params['MASTER_STILTS.PY']['match_tol'])
    modo_ms   = str(params['MASTER_STILTS.PY']['modo_ms'])

    return match_tol, modo_ms

def get_master_match():
    min_epochs_mm = int(params['MASTER_MATCH_ID.PY']['min_epochs_mm'])
    nprocs_mmi    = int(params['MASTER_MATCH_ID.PY']['nprocs_mmi'])
    return min_epochs_mm, nprocs_mmi

def get_mastercat():
    min_epochs, min_mag, max_mag, max_err, iteraciones, iteracion2, nrefstars_mc, nprocs_mc = params['MASTERCAT.PY'].values()
    min_epochs, min_mag, max_mag, max_err, iteraciones, nrefstars_mc, nprocs_mc = np.array([min_epochs, min_mag, max_mag, max_err, iteraciones, nrefstars_mc, nprocs_mc]).astype(float)

    iteracion2 = str(iteracion2)
    return min_epochs, min_mag, max_mag, max_err, iteraciones, iteracion2, nrefstars_mc, nprocs_mc

def get_match1a1():
    modo_ma = str(params['MATCH_1a1.PY']['modo_ma'])
    tol     = float(params['MATCH_1a1.PY']['tol'])
    nprocs_m1a1 = int(params['MATCH_1a1.PY']['nprocs_m1a1'])

    return modo_ma, tol, nprocs_m1a1

def get_pm1a1():
    nframes, nbins, limplotpm, nprocs_pm, sig_iter, nsigma, weight = params['PM_1a1.PY'].values()
    nframes, nbins, limplotpm, nprocs_pm, sig_iter, nsigma = np.array([nframes, nbins, limplotpm, nprocs_pm, sig_iter, nsigma]).astype(float)
    nprocs_pm = int(nprocs_pm)
    sig_iter  = int(sig_iter)
    nsigma    = int(nsigma)

    weight = str(weight) == 'True'

    return nframes, nbins, limplotpm, nprocs_pm, sig_iter, nsigma, weight

def get_XYtoRADEC():
    nprocs = np.array(params['XYtoRADEC.PY'].values()).astype(float)

    return nprocs

#ProgressBar con multiprocess (elegir numero de procesadores)
from astropy.utils.console import ProgressBar
import multiprocessing
def barra(funcion, items, cpus):
    results = []
    with ProgressBar(len(items)) as bar:
        p  = multiprocessing.Pool(processes=cpus)

        for jj, result in enumerate(p.imap(funcion, items)):
            bar.update(jj)
            results.append(result)
        p.close()
        p.join()
    return results
