#Crea un master frame usando como referencia(s) epocas dadas por el usuario

import numpy as np
import glob
import sys
from astropy import wcs
from astropy.utils.console import ProgressBar, color_print
from astropy.io import fits
from scipy.optimize import curve_fit

#Para ignorar los warnings de WCS
import warnings
warnings.filterwarnings("ignore")
