import pluto
import tools
import numpy as np
import matplotlib.pyplot as plt
from celluloid import Camera
import os 

import UI_helper
import Navigation_helper
import Data_parser_helper
from Global_variables import *
import plot_params
import plotter_helper

def plot_one_disk_accretion() -> None:
    data_name = UI_helper.selectDataToPlot()
    plot_disk_accretion(data_name)

def plot_one_planet_accretion() -> None:
    data_name = UI_helper.selectDataToPlot()
    plot_planet_accretion(data_name)

def plot_disk_accretion(data_name: str) -> None:
    plot_params.square()
    directories = Navigation_helper.Directories(data_name)
    
    # obj, obj_des = UI_helper.selectObjectToPlot(directories.out_dir)
    (time, _, _, _, _, _, _, _, _, accretion_list) = Data_parser_helper.get_averages_data(data_name) 
    total_accretion = [
    sum([accretion_list[j][i] for j in range(len(accretion_list))]) 
    for i in range(len(accretion_list[0]))
    ]

    fig = plt.figure()
    plt.plot(time, total_accretion)
    fig.tight_layout()

    f_name = '{}disk_accretion'.format(directories.plots_dir)
    save_path = plotter_helper.define_save_plot(f_name)
    fig.savefig(save_path)
    plt.close(fig)

def plot_planet_accretion(data_name: str) -> None:
    plot_params.square()
    directories = Navigation_helper.Directories(data_name)
    object = UI_helper.selectObjectToPlot(directories.out_dir)
    
    (
        time,
        a,
        e,
        period,
        mass,
    ) = Data_parser_helper.getNbodyInformation_dat(data_name, object.id) 
    if not isinstance(mass, np.ndarray): raise Exception('no mass data available')
    delta_mass = []
    fig = plt.figure()
    plt.plot(time, mass)
    plt.xlabel('Time (binary orbits)')
    plt.ylabel('Planet mass')
    plt.legend()
    plt.grid()
    fig.tight_layout()

    f_name = '{}obj{}_accretion'.format(directories.plots_dir, object.id)
    save_path = plotter_helper.define_save_plot(f_name, 'png')
    fig.savefig(save_path)
    plt.close(fig)

if __name__ == '__main__':
    plotters = [plot_one_disk_accretion, plot_one_planet_accretion]
    func_index = UI_helper.selectFunctionsToRun(plotters)
    eval('{}()'.format(plotters[func_index].__name__))