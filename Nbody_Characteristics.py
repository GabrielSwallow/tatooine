# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 15:42:15 2022

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

a_bin = 1 #units of binary separation
size = 7.5 #limits of graph
Rmax = 70 #unused

cart = False # 'cart grid'
logsc = False # log scale
nbody = True # nbody integrater used
nts   = 100 # output time step
var   = 'rho' # 'rho, vx1, vx2, prs'

name = 'Kep-47' #name of the simulation
fs = 14 #font size

def plot_one() -> None:
    data_name = UI_helper.selectDataToPlot()
    plot(data_name)

def plot(data_name: str):
    obj = UI_helper.selectObjectToPlot()
    
    if obj ==1:
        obj_des = 'B'
    elif obj ==2:
        obj_des = 'b'
    elif obj ==3:
        obj_des = 'c'
    elif obj ==4:
        obj_des = 'd'
        
    all_data_dir = '../Data/'
    
    data_parent_dir = all_data_dir + data_name
    out_dir = all_data_dir + data_name + '/out'
    
    fig, axs = plt.subplots(4,1, sharex= 'all')
    data_parent_dir = all_data_dir + data_name
    out_dir = all_data_dir + data_name + '/out'
    plots_dir = data_parent_dir + '/Plots/'
    nbody_elements = out_dir + '/nbody_orbital_elements.out'
    plutolog = out_dir + '/pluto.log'
    
    dbl_txt = open(plutolog)
    dbl_strs = dbl_txt.readlines()
    dbl_str = dbl_strs[79]
    dbl_list = list(map(int, re.findall('\d+', dbl_str)))
    dbl_time = dbl_list[3]
    dbl_txt.close()
    
    (unfiltered_time,
    object_id,
    unfiltered_a,
    unfiltered_e,
    unfiltered_omega,
    unfiltered_anomaly) = np.loadtxt(nbody_elements,
                                     usecols=(0,1,2,3,6,7),
                                     unpack = True)
    time = unfiltered_time[object_id == obj]
    time *= dbl_time
    a = unfiltered_a[object_id == obj]
    e = unfiltered_e[object_id == obj]
    omega = unfiltered_omega[object_id == obj]
    anomaly = unfiltered_anomaly[object_id == obj]
    
    axs[0].plot(time, a)
    axs[0].set_title('Semi-Major Axis [$\mathrm{a_{bin}}$]')
    
    axs[1].plot(time, e)
    axs[1].set_title('Eccentricity')
    
    axs[2].plot(time, omega)
    axs[2].set_title('Angular Distance of Pericenter [$\mathrm{Rad}$]')
    
    axs[3].plot(time, anomaly)
    axs[3].set_title('True Anomaly [$\mathrm{Rad}$]')
    axs[3].set(xlabel = 'Time [$\mathrm{T_{bin}}$]')
    fig.suptitle('Orbital Elements of Kepler-47{}'.format(obj_des))
    
    fig.tight_layout()
    
    save_path = '{}{}_obj{}_orbital_elements.png'.format(plots_dir, data_name, obj)
    
    print('Saving plot in {0}'.format(save_path))
    fig.savefig(save_path)
    plt.close(fig)

if __name__ == '__main__':
    plot_one()
