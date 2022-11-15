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
from Global_variables import *


def plot_one() -> None:
    data_name = UI_helper.selectDataToPlot()
    plot(data_name)

def plot(data_name: str) -> None:
    obj = UI_helper.selectObjectToPlot()
    data_parent_dir = all_data_dir + data_name
    out_dir = all_data_dir + data_name + '/out'
    plots_dir = data_parent_dir + '/Plots/'

    nbody_elements_filename = out_dir + '/nbody_orbital_elements.out'
    pluto_log_filename = out_dir + '/pluto.log'
    
    if obj == 0:
        print('\nCannot select obj = 0 for Nbody characteristics. Try again.')
        # raise Exception('Cannot select obj = 0 for Nbody characteristics')
        plot(data_name)
        return
    if obj == 1:
        obj_des = 'B'
    elif obj == 2:
        obj_des = 'b'
    elif obj == 3:
        obj_des = 'c'
    elif obj == 4:
        obj_des = 'd' 
    
    

    dbl_txt = open(pluto_log_filename)
    dbl_strs = dbl_txt.readlines()
    dbl_str = dbl_strs[79]
    dbl_list = list(map(int, re.findall('\d+', dbl_str)))
    dbl_time = dbl_list[3]
    dbl_txt.close()
    
    (
        unfiltered_time,
        object_id,
        unfiltered_a,
        unfiltered_e,
        unfiltered_omega,
        unfiltered_anomaly,
    ) = np.loadtxt(
        nbody_elements_filename, usecols=(0,1,2,3,6,7), unpack = True
        )

    time = unfiltered_time[object_id == obj]
    if len(time) == 0: 
        print('\nMax object_id = {}. You selected {}. Try again'.format(int(max(object_id)), obj))
        plot(data_name)
        return
    time *= dbl_time
    a = unfiltered_a[object_id == obj]
    e = unfiltered_e[object_id == obj]
    omega = unfiltered_omega[object_id == obj]
    anomaly = unfiltered_anomaly[object_id == obj]
    
    fig, axs = plt.subplots(4,1, sharex= 'all')

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
    repeated_plots = 0
    while(os.path.isfile(save_path)):
        save_path = '{}{}_obj{}_orbital_elements({}).png'.format(plots_dir, data_name, obj, repeated_plots)
        repeated_plots += 1
    
    print('Saving plot in {0}'.format(save_path))
    fig.savefig(save_path)
    plt.close(fig)

if __name__ == '__main__':
    plot_one()