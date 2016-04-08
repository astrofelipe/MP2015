from configparser import ConfigParser
import numpy as np

params = ConfigParser()
params.read('zparams_pm.py')

def get_script():
    radio, itera = np.array(params['SCRIPT.PY'].values()).astype(int)
    output  = str(params['TLINEAL_1a1.PY']['output'])
    refer   = str(params['TLINEAL_1a1.PY']['refer'])
    nframes = int(params['PM_1a1.PY']['nframes'])
    min_ep  = int(params['VPDHmag.PY']['min_ep'])

    return radio, itera, output, refer, nframes, min_ep

def get_VPDHmag():
    limplotvp, magl, magh, delta, min_ep, min_nei, sigma_err = np.array(params['VPDHmag.PY'].values()).astype('float')
    mags = magl, magh

    return limplotvp, mags, delta, min_ep, min_nei, sigma_err

def get_tlineal():
    nrefstars_tl, min_nei, rad_int, rad_ext, output, refer, sort_mag, \
    local, ma1, ma2, mr1, mr2, mp1, mp2, rad_ref, x0, y0, lim, plot_ep, plot_del_ep, plot_del_xy = params['TLINEAL_1a1.PY'].values()

    nrefstars_tl, min_nei, rad_int, rad_ext, ma1, ma2, mr1, mr2, mp1, mp2, rad_ref, x0, y0, lim = \
    np.array([nrefstars_tl, min_nei, rad_int, rad_ext, ma1, ma2, mr1, mr2, mp1, mp2, rad_ref, x0, y0, lim]).astype(float)

    refer  = str(refer)
    output = str(output)

    sort_mag    = str(sort_mag) == 'True'
    local       = str(local) == 'True'
    plot_ep     = str(plot_ep) == 'True'
    plot_del_ep = str(plot_del_ep) == 'True'
    plot_del_xy = str(plot_del_xy) == 'True'

    return nrefstars_tl, min_nei, rad_int, rad_ext, output, refer, sort_mag, \
    local, ma1, ma2, mr1, mr2, mp1, mp2, rad_ref, x0, y0, lim, plot_ep, plot_del_ep, plot_del_xy

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
    min_epochs, min_mag, max_mag, max_err, iteraciones, iteracion2, nrefstars_mc = params['MASTERCAT.PY'].values()
    min_epochs, min_mag, max_mag, max_err, iteraciones, nrefstars_mc = np.array([min_epochs, min_mag, max_mag, max_err, iteraciones, nrefstars_mc]).astype(float)

    iteracion2 = str(iteracion2)
    return min_epochs, min_mag, max_mag, max_err, iteraciones, iteracion2, nrefstars_mc

def get_match1a1():
    modo_ma = str(params['MATCH_1a1.PY']['modo_ma'])
    tol     = float(params['MATCH_1a1.PY']['tol'])

    return modo_ma, tol

def get_pm1a1():
    nframes, nbins, limplotpm = np.array(params['PM_1a1.PY'].values()).astype(float)

    return nframes, nbins, limplotpm

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
