# -*- coding: utf-8 -*-

import numpy as np

###funciones

def funcion1(x):
    return x**2-2

def funcion2(x):
    return (x**5)-6.6*(x**4)+5.12*(x**3)+21.312*(x**2)-38.016*(x)+17.28

def funcion3(x):
    return (x-1.5)*(np.e**(-4*((x-1.5)**2)))

def funcion1p(x):
    return x*2

def funcion2p(x):
    return 5*(x**4)-26.4*(x**3)+15.36*(x**2)+42.624*(x)-38.016

def funcion3p(x):
    return (x-1.5)*(np.e**(-4*((x-1.5)**2)))*(-8*x+12.0)+np.e**((-4*(x-1.5))**2) 

def funcion1pp(x):
    return 2

def funcion2pp(x):
    return 20*(x**3)-79.2*(x**2)+30.72*x+42.624

def funcion3pp(x):
    return (-24*x+(x-1.5)*((8*x-12)**2)+36)*np.e**(-4*((x-1.5)**2))