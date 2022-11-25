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


def plot_one_using_dat() -> None:
    data_name = UI_helper.selectDataToPlot()
    plot_using_dat(data_name)

def plot_one_using_out() -> None:
    data_name = UI_helper.selectDataToPlot()
    plot_using_out(data_name)

def plot_using_dat(data_name: str) -> None:
    directories = Navigation_helper.Directories(data_name)
    
    obj, obj_des = UI_helper.selectObjectToPlot(directories.out_dir)
    
    (
        time,
        a,
        e,
        period,
    ) = Data_parser_helper.getNbodyInformation_dat(directories.out_dir, obj) 

    fig, axs = plt.subplots(2,1, sharex= 'all')

    fin = 100
    axs[0].plot(time[0:fin], a[0:fin])
    axs[0].set_title('Semi-Major Axis [$\mathrm{a_{bin}}$]')
    
    axs[1].plot(time[0:fin], e[0:fin])
    axs[1].set_title('Eccentricity')
    
    # axs[2].plot(time, period)
    # axs[2].set_title('period [$\mathrm{Rad}$]')
    axs[1].set(xlabel = 'Time [$\mathrm{T_{bin}}$]')
    fig.suptitle('Orbital Elements of Kepler-47{}'.format(obj_des))
    
    fig.tight_layout()
    
    save_path = '{}{}_obj{}_orbital_elements_dat.png'.format(directories.plots_dir, data_name, obj)
    repeated_plots = 0
    while(os.path.isfile(save_path)):
        save_path = '{}{}_obj{}_orbital_elements_dat({}).png'.format(directories.plots_dir, data_name, obj, repeated_plots)
        repeated_plots += 1
    
    print('Saving plot in {0}'.format(save_path))
    fig.savefig(save_path)
    plt.close(fig)

def plot_using_out(data_name: str) -> None:
    directories = Navigation_helper.Directories(data_name)

    obj, obj_des = UI_helper.selectObjectToPlot(directories.out_dir)
    
    (
        time,
        a,
        e,
        omega,
        anomaly,
    ) = Data_parser_helper.getNbodyInformation_out(directories.out_dir, obj) 
    
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
    
    save_path = '{}{}_obj{}_orbital_elements.png'.format(directories.plots_dir, data_name, obj)
    repeated_plots = 0
    while(os.path.isfile(save_path)):
        save_path = '{}{}_obj{}_orbital_elements({}).png'.format(directories.plots_dir, data_name, obj, repeated_plots)
        repeated_plots += 1
    
    print('Saving plot in {0}'.format(save_path))
    fig.savefig(save_path)
    plt.close(fig)

if __name__ == '__main__':
    plotters = [plot_one_using_dat, plot_one_using_out]
    func_index = UI_helper.selectFunctionsToRun(plotters)
    eval('{}()'.format(plotters[func_index].__name__))