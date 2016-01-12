#Los valores de string deben ser entregados sin comillas, salvo casos como "RA DEC" cuando
#se quiere hacer un match (stilts) con mas de una columna

[CMD.PY]
#Modo para seleccionar los limites (manual o auto)
cmd_modo = manual

#Limites en color
col1     = -2
col2     = 2

#Limite en magnitud
mag1     = 20
mag2     = 10

[MASTER_STILTS.PY]
#Tolerancia para el modo a usar (ID no usa este parametro)
match_tol     = .3

#Columnas para hacer el match (ej "RA DEC" o ID)
modo_ms       = ID

[MASTER_MATCH_ID.PY]
#Minimo de epocas para considerar una estrella
min_epochs_mm = 2

[MASTERCAT.PY]
#Numero minimo de epocas para considerar la estrella (mastercat, PM_1a1, VPDHmag y master_match_id)
min_epochs  = 2

#Magnitud minima y maxima para estrellas a transformar
min_mag     = 11
max_mag     = 14

#Error maximo a considerar
max_err     = .05

#Numero de iteraciones a realizar
iteraciones = 3

#Como proceder para la segunda iteracion en adelante (global o local)
iteracion2   = global

#Numero de vecinos locales (contando a la estrella propia) si iteracion2 = local
nrefstars_mc = 51

[MATCH_1a1.PY]
#ID o "RA DEC" (comillas para varias columnas)
modo_ma   = ID

#Tolerancia match (caso RA DEC)
tol       = 0.3

[TLINEAL_1a1.PY]
#Numero de refstars(+1) deseadas
nrefstars_tl = 51

#Minimo de vecinos para considerar (sino tira 888)
min_nei   = 4

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
ma2 = 14.0

#Corte en magnitud para las refstars
mr1 = 11.0
mr2 = 14.0

#Corte en magnitud para plotear las refstars
mp1 = 11.0
mp2 = 12.0

#Radio (pix) dentro del cual una refstar se considera para plot
rad_ref = 99999999

#Coordenadas centrales del circulo a considerar
x0 = 1352
y0 = 554

#Limites del plot en pixeles (cuadrado, por eso es uno)
lim = 2

#Plot delta vs epocas y delta vs coor_x o coor_y
plot_del_ep = True
plot_del_xy = True

[PM_1a1.PY]
#Numero minimo de epocas en que debe estar la estrella
nframes   = 2

#Delta del bin
nbins     = 1

#Limites del plot VPD en mas/yr
limplotpm = 30

[VPDHmag.PY]
#Limites del plot (bineado) en mas/yr
limplotvp = 30

#Limite de magnitudes para plotear
mags_l  = 11
mags_h  = 14

#Intervalo de magnitudes
delta   = .5

#Minimo de epocas en que debe estar una estrella para considerar
min_ep  = 2

[SCRIPT.PY]
#Radio en arcsec para seleccionar las nuevas refstars de PM_final.dat
radio = 3

#Numero de iteraciones
itera = 3
