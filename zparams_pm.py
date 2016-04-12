#Los valores de string deben ser entregados sin comillas, salvo casos como "RA DEC" cuando
#se quiere hacer un match (stilts) con mas de una columna
[XYtoRADEC.PY]
nprocs = 1


[CMD.PY]
#Modo para seleccionar los limites (manual o auto)
cmd_modo = manual

#Columnas para usar en match (ID o RA DEC)
match = ID

#Limites en color
col1     = -2
col2     = 2

#Limite en magnitud
mag1     = 20
mag2     = 10

#Genera PDF
cmd_pdf  = True

[MASTER_STILTS.PY]
#Tolerancia para el modo a usar (ID no usa este parametro)
match_tol     = .3

#Columnas para hacer el match (ej RA DEC o ID)
modo_ms       = ID

[MASTER_MATCH_ID.PY]
#Minimo de epocas para considerar una estrella
min_epochs_mm = 8
nprocs_mmi    = 1

[MASTERCAT.PY]
#Numero minimo de epocas para considerar la estrella (mastercat, PM_1a1, VPDHmag y master_match_id)
min_epochs  = 8

#Magnitud minima y maxima para estrellas a transformar
min_mag     = 11
max_mag     = 15

#Error maximo a considerar
max_err     = .05

#Numero de iteraciones a realizar
iteraciones = 3

#Como proceder para la segunda iteracion en adelante (global o local)
iteracion2   = global

#Numero de vecinos locales (contando a la estrella propia) si iteracion2 = local
nrefstars_mc = 51

nprocs_mc   = 1


[MATCH_1a1.PY]
#ID o "RA DEC" (comillas para varias columnas)
modo_ma   = ID

#Tolerancia match (caso RA DEC)
tol       = 0.3

#Numero procesadores
nprocs_m1a1 = 1

[TLINEAL_1a1.PY]
#Numero de refstars(+1) deseadas
nrefstars_tl = 51

#Minimo de vecinos para considerar (sino tira 888)
min_nei   = 5

#Radio interior y exterior. Exterior = 0 -> selecciona nrefstars mas cercanas
rad_int   = 1
rad_ext   = 300

#PDF de Output (reemplazara sin avisar!)
output    = test

#Catalogo con las estrellas de referencia
refer     = refstars.gc

#Sort refstars segun mag y toma las nrefstars mas brillantes
sort_mag  = False

#True para usar transformaciones locales, False para usar global
local     = True

#Corte en magnitud para considerar estrellas a analizar
ma1 = 11.0
ma2 = 19.0

#Corte en magnitud para las refstars
mr1 = 11.0
mr2 = 15.0

#Corte en magnitud para plotear las refstars
mp1 = 11.0
mp2 = 14.0

#Radio (pix) dentro del cual una refstar se considera para plot
rad_ref = 99999999

#Coordenadas centrales del circulo a considerar
x0 = 1352
y0 = 554

#Limites del plot en pixeles (cuadrado, por eso es uno)
lim = 2

#Plot de epocas
plot_ep = False

#Plot delta vs epocas y delta vs coor_x o coor_y
plot_del_ep = False
plot_del_xy = False

#Numero de procesadores
nprocs_tl = 1

[PM_1a1.PY]
#Numero minimo de epocas en que debe estar la estrella
nframes   = 8

#Delta del bin
nbins     = 1

#Limites del plot VPD en mas/yr
limplotpm = 20

#Numero procesadores
nprocs_pm = 1

[VPDHmag.PY]
#Limites del plot (bineado) en mas/yr
limplotvp = 20

#Limite de magnitudes para plotear
mags_l  = 9
mags_h  = 19

#Intervalo de magnitudes
delta   = 1

#Minimo de epocas en que debe estar una estrella para considerar
min_ep  = 8

#Minimo de vecinos con los que hizo la transformacion
min_nei = 5

#Numero de sigmas para hacer rejection por errores (modulo pmxe y pmye)
sigma_err = 3

[SCRIPT.PY]
#Radio en arcsec para seleccionar las nuevas refstars de PM_final.dat
radio = 3

#Numero de iteraciones
itera = 3
