# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 17:30:13 2022

@author: epick
"""

import pluto
import tools
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib as mpl
import matplotlib.colors as colors
import argparse as argp
import os
import UI_helper
import Navigation_helper
import re
import scipy.optimize as opt

def ellipse(phi, a, e, phi0):
    return a*(1-e**2)/(1+e*np.cos(phi - phi0))

def x_coord(r, phi):
    return r*np.cos(phi)

def y_coord(r, phi):
    return r*np.sin(phi)

def r_coord(x, y):
    return np.sqrt(x**2 + y**2)

def ellipse_find(R, Phi, Sigma):
    dimensions = np.shape(R)
    phi_range = dimensions[0]
    r_range = dimensions[1]
    
    radii = []
    phis = []
    
    for i in range(0, phi_range - 1):
        phi = Phi[i,0]
        phis.append(phi)
        sigma = Sigma[i,:]
        r = R[i,:]
        sigma_max = np.amax(sigma)
        for j in range(0, r_range - 1):
            if sigma[j] < 0.1 * sigma_max:
                continue
            else:
                break
        radii.append(r[j])
    
    radii = np.array(radii)
    phis = np.array(phis)
    x = x_coord(radii, phis)
    y = y_coord(radii, phis)
    
    x0 = np.mean(x)
    y0 = np.mean(y)
    
    x_fix = x - x0
    y_fix = y - y0
    
    r_fix = r_coord(x_fix, y_fix)
    
    popt, _ = opt.curve_fit(ellipse, phis, r_fix, p0= [5.0, 0.3, 0.1])
    
    a = popt[0]
    e = popt[1]
    omega = popt[2]
    
    return a, e, omega, x0, y0
