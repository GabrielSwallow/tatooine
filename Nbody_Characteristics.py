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
import Data_parser_helper
import plot_params


def plot_one_using_dat() -> None:
    data_name = UI_helper.selectDataToPlot()
    plot_using_dat(data_name)

def plot_one_using_out() -> None:
    data_name = UI_helper.selectDataToPlot()
    plot_using_out(data_name)

def plot_using_dat(data_name: str) -> None:
    plot_params.square()
    directories = Navigation_helper.Directories(data_name)
    
    obj, obj_des = UI_helper.selectObjectToPlot(directories.out_dir)
    n_min, n_max = UI_helper.selectPlottingRange(directories.out_dir)
    t_min = n_min * nts
    t_max = n_max * nts
    
    (
        time,
        a,
        e,
        period,
        mass,
    ) = Data_parser_helper.getNbodyInformation_dat(data_name, obj) 

    i_min, i_max = tools.time_split(time, t_min, t_max)

    fig, axs = plt.subplots(2,1, sharex= 'all')

    fin = len(time) - 1
    axs[0].plot(time[i_min:i_max], a[i_min:i_max])
    axs[0].set_title('Semi-Major Axis [$\mathrm{a_{bin}}$]')
    
    axs[1].plot(time[i_min:i_max], e[i_min:i_max])
    axs[1].set_title('Eccentricity [$\mathrm{e}$]')
    
    # axs[2].plot(time, period)
    # axs[2].set_title('period [$\mathrm{Rad}$]')
    axs[1].set(xlabel = 'Time [$\mathrm{T_{bin}}$]')
    fig.suptitle('Orbital Elements of Kepler-47{}'.format(obj_des))
    
    fig.tight_layout()
    
    save_path = '{}{}_obj{}_orbital_elements_dat_{}-{}.png'.format(directories.plots_dir, data_name, obj, n_min, n_max)
    repeated_plots = 0
    while(os.path.isfile(save_path)):
        save_path = '{}{}_obj{}_orbital_elements_dat_{}-{}({}).png'.format(directories.plots_dir, data_name, obj, n_min, n_max, repeated_plots)
        repeated_plots += 1
    
    print('Saving plot in {0}'.format(save_path))
    fig.savefig(save_path)
    plt.close(fig)

def plot_using_out(data_name: str) -> None:
    plot_params.square()
    directories = Navigation_helper.Directories(data_name)

    obj, obj_des = UI_helper.selectObjectToPlot(directories.out_dir)
    n_min, n_max = UI_helper.selectPlottingRange(directories.out_dir)
    
    (
        time,
        a,
        e,
        omega,
        anomaly,
    ) = Data_parser_helper.getNbodyInformation_out(directories.data_name, obj) 
    
    fig, axs = plt.subplots(4,1, sharex= 'all')

    axs[0].plot(time[n_min:n_max], a[n_min:n_max])
    axs[0].set_title('Semi-Major Axis [$\mathrm{a_{bin}}$]')
    
    axs[1].plot(time[n_min:n_max], e[n_min:n_max])
    axs[1].set_title('Eccentricity')
    
    axs[2].plot(time[n_min:n_max], omega[n_min:n_max])
    axs[2].set_title('Angular Distance of Pericenter [$\mathrm{Rad}$]')
    
    axs[3].plot(time[n_min:n_max], anomaly[n_min:n_max])
    axs[3].set_title('True Anomaly [$\mathrm{Rad}$]')
    axs[3].set(xlabel = 'Time [$\mathrm{T_{bin}}$]')
    fig.suptitle('Orbital Elements of Kepler-47{}'.format(obj_des))
    
    fig.tight_layout()
    
    save_path = '{}{}_obj{}_orbital_elements_{}-{}.png'.format(directories.plots_dir, data_name, obj, n_min, n_max)
    repeated_plots = 0
    while(os.path.isfile(save_path)):
        save_path = '{}{}_obj{}_orbital_elements_{}-{}({}).png'.format(directories.plots_dir, data_name, obj, n_min, n_max, repeated_plots)
        repeated_plots += 1
    
    print('Saving plot in {0}'.format(save_path))
    fig.savefig(save_path)
    plt.close(fig)

if __name__ == '__main__':
    plotters = [plot_one_using_dat, plot_one_using_out]
    func_index = UI_helper.selectFunctionsToRun(plotters)
    eval('{}()'.format(plotters[func_index].__name__))