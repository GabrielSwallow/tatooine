import os
from typing import Tuple
import Navigation_helper
import Data_parser_helper
from Global_variables import *
import plot_params
import matplotlib.pyplot as plt
import UI_helper
from collections.abc import Callable
import numpy as np

# from menu import Menu
# import simple_term_menu as stm

def plot_multiple_data_sets_overlayed(
    many_data_to_plot: list[data_id], 
    plotter: Callable[[plt.Axes, str, str, list], None], 
    plotter_args: list = None
    ) -> None:
    plot_params.square()

    fig, ax = plt.subplots(1,1)
    for data_to_plot in many_data_to_plot:
        if plotter_args == None:
            plotter(ax, data_to_plot.name, data_to_plot.legend_name)
        else:
            plotter(ax, data_to_plot.name, data_to_plot.legend_name, plotter_args)
    
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


def plot_instability_zone(ax: plt.Axes) -> None:
    radius = 2.31
    theta = np.arange(0, 2*np.pi, 0.01)
    ax.plot(radius*np.cos(theta), radius*np.sin(theta), '--w', lw=1.5)

def plot_Kep47b_for_line_plot(ax: plt.Axes, t_min: float, t_max: float) -> None:
    K_47b_radius = objects_in_kep47[2].radius
    ax.plot([t_min,t_max], [K_47b_radius,K_47b_radius], 'g--', linewidth = 2, label='kep47b final orbit')

def plot_Kepler_47_planets_for_twoD_sigma(ax: plt.Axes) -> None:
    K_47b_radius = objects_in_kep47[2].radius
    K_47c_radius = objects_in_kep47[3].radius
    K_47d_radius = objects_in_kep47[4].radius
    theta = np.arange(0, 2*np.pi, 0.01)
    ax.plot(K_47b_radius*np.cos(theta), K_47b_radius*np.sin(theta), '#fffdcfff', lw=3)
    ax.plot(K_47d_radius*np.cos(theta), K_47d_radius*np.sin(theta), '#fffdcfff', lw=3)
    ax.plot(K_47c_radius*np.cos(theta), K_47c_radius*np.sin(theta), '#fffdcfff', lw=3)