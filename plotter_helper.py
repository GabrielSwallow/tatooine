import os
from typing import Tuple
import Navigation_helper
import Data_parser_helper
from Global_variables import *
import plot_params
import matplotlib.pyplot as plt
import UI_helper
from collections.abc import Callable

# from menu import Menu
# import simple_term_menu as stm

def plot_multiple_data_sets_overlayed(
    many_data_to_plot: list[data_id], 
    plotter: Callable[[plt.Figure, str, str, list], None], 
    plotter_args: list = None
    ) -> None:
    plot_params.square()

    fig = plt.figure()
    for data_to_plot in many_data_to_plot:
        if plotter_args == None:
            plotter(fig, data_to_plot.name, data_to_plot.legend_name)
        else:
            plotter(fig, data_to_plot.name, data_to_plot.legend_name, plotter_args)
    
    plot_name = UI_helper.name_the_plot() + '_'

    save_path = '{}disc_eccentricity_averages_{}.png'.format(global_plots_dir+plot_name, var)
    repeated_plots = 1
    while(os.path.isfile(save_path)):
        save_path = '{}disc_eccentricity_averages_{}({}).png'.format(global_plots_dir+plot_name, var, repeated_plots)
        repeated_plots += 1
    print('Saving plot in {0}'.format(save_path))
    plt.legend()
    plt.xlabel('Time [$\mathrm{T_{bin}}$]')
    plt.savefig(save_path)
    plt.close()

def plot_multiple_data_sets_overlayed_subplots(
    many_data_to_plot: list[data_id], 
    plotter: Callable[[plt.Figure, plt.Axes, str, str, list], None], 
    plotter_args: list = None
    ) -> None:
    plot_params.square()

    fig, axs = plt.subplots(2,1, sharex = 'all')

    for data_to_plot in many_data_to_plot:
        if plotter_args == None:
            plotter(fig, axs, data_to_plot.name, data_to_plot.legend_name)
        else:
            plotter(fig, axs, data_to_plot.name, data_to_plot.legend_name, plotter_args)

    
    plot_name = UI_helper.name_the_plot() + '_'

    save_path = '{}disc_eccentricity_averages_{}.png'.format(global_plots_dir+plot_name, var)
    repeated_plots = 1
    while(os.path.isfile(save_path)):
        save_path = '{}disc_eccentricity_averages_{}({}).png'.format(global_plots_dir+plot_name, var, repeated_plots)
        repeated_plots += 1
    print('Saving plot in {0}'.format(save_path))
    plt.legend()
    plt.xlabel('Time [$\mathrm{T_{bin}}$]')
    plt.savefig(save_path)
    plt.close()