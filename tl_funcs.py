from __future__ import division, print_function
import numpy as np

def linear(coords, a, b, c):
    x, y = coords
    return a*x + b*y + c

def quadratic(coords, a, b, c, d, e ,f):
    x, y = coords
    return a*x + b*y + c*x**2 + d*y**2 + e*x*y + f
