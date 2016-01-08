#Archivo de parametros para las otras rutinas

#MASTER_STILTS.PY
match_tol     = .3           #Tolerancia para el modo a usar (ID no usa este parametro)
modo_ms       = '"RA DEC"'    #Columnas para hacer el match (ej 'RA DEC' o 'ID')
#modo          = 'ID'         #Columnas para hacer el match (ej 'RA DEC' o 'ID')

#MASTER_MATCH_ID.PY
min_epochs_mm = 5

#MASTERCAT.PY
min_epochs  = 2  #Numero minimo de epocas para considerar la estrella (mastercat, PM_1a1, VPDHmag y master_match_id)
min_mag     = 11     #Magnitud minima para estrellas a transformar
max_mag     = 14     #Magnitud maxima...
max_err     = .05    #Error maximo a considerar
iteraciones = 5
iteracion2  = 'global'   #global o local

nrefstars_mc   = 51 #Numero de vecinos locales (contando a la estrella propia) si iteracion2=='local'

#MATCH_1a1.PY
modo_ma   = '"RA DEC"'            #'ID' o '"RA DEC"' (comillas para varias columnas)
tol       = 0.3             #Tolerancia match (caso RA DEC)

#TLINEAL_1a1.PY
nrefstars_tl = 51         #numero de refstars(+1) deseadas
min_nei   = 4         #Minimo de vecinos para considerar (sino tira 888)
rad_int   = 1         #Radio interior
rad_ext   = 300         #Radio exterior (0 -> selecciona nrefstars mas cercanas)
output    = 'test'     #PDF de Output (reemplazara sin avisar!)
refer     = 'refstars.gc' #Catalogo con las estrellas de referencia
sort_mag  = False      #Sort refstars segun mag y toma las nrefstars mas brillantes
local     = True           #True para usar transf locales, False para usar global
ma1,ma2   = 11.0,14.0      #Corte en magnitud para considerar estrellas a analizar
mr1,mr2   = 11.0,14.0      #Corte en magnitud para las refstars
mp1,mp2   = 11.0,12.0     #Corte en magnitud para plotear las refstars

rad_ref  = 99999999       #Radio (pix) dentro del cual una refstar se considera para plot
x0,y0    = 1352,554       #Coordenadas centrales del circulo a considerar

lim    = 2         #Limites del plot en pixeles (cuadrado, por eso es uno)

plot_del_ep = True     #Plot delta vs epocas
plot_del_xy = True      #Plot delta vs coor_x o coor_y

#PM_1a1.PY
nframes   = 2    #Numero minimo de epocas en que debe estar la estrella
nbins     = 1    #Delta del bin
limplotpm = 30   #Limites del plot VPD en mas/yr

#VPDHmag.PY
limplotvp = 30          #Limites del plot (bineado) en mas/yr
mags      = 11, 14      #Limite de magnitudes para plotear
delta     = .5          #Intervalo de magnitudes
min_ep    = 2           #Minimo de epocas en que debe estar una estrella para considerar

#SCRIPT.PY
radio = 3   #Radio en arcsec para seleccionar las nuevas refstars de PM_final.dat
itera = 3   #Numero de iteraciones
