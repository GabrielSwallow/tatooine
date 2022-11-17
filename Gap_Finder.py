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

#%%

def ellipse(phang, a, e):
    return a*(1-e**2)/(1+e*np.cos(phang))

var = 'rho'

n = 19

datafile = '../Data/141122_0/out'

data = pluto.Pluto(datafile)

R, Phi = np.meshgrid(data.grid['faces'][0], data.grid['faces'][1])
Sigma = data.primitive_variable(var, n)[0,:,:]

r = R[1,:]

dimensions = np.shape(R)

phirange = dimensions[0]

radii = []
angles = []

for i in range(0,phirange-1):
    phi = Phi[i,0]
    sigma = Sigma[i,:]
    r = R[i,:]
    smax = np.amax(sigma)
    index = np.where(sigma==smax)
    radii.append(r[index[0]])
    angles.append(phi)
    
radii_ideal = ellipse(angles, 5.18, 0.39)

#%%

x=[]
y=[]
xid = []
yid = []
for i in range(0,phirange-1):
    x.append(radii[i]*np.cos(angles[i]))
    y.append(radii[i]*np.sin(angles[i]))
    xid.append(radii_ideal[i]*np.cos(angles[i]))
    yid.append(radii_ideal[i]*np.sin(angles[i]))

#%%

fig = plt.figure()
ax = fig.add_subplot()
plt.plot(x,y)
plt.plot(xid,yid)
plt.xlim([-8,8])
plt.ylim([-8,8])
ax.set_aspect('equal', adjustable='box')
plt.show()

#%%

radii = np.array(radii)
angles = np.array(angles)

popt, _ = opt.curve_fit(ellipse, angles, radii, p0 = np.array([5.18,0.39]))
