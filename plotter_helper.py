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
from tools import Unit_conv

# from menu import Menu
# import simple_term_menu as stm

def define_save_plot(plots_dir: str, file_name: str, extension: str = 'png', overwrite: bool = overwrite_plots):
    save_path = '{}{}.{}'.format(plots_dir, file_name, extension)
    if not overwrite_plots:
        repeated_plots = 0
        while(os.path.isfile(save_path)):
            save_path = '{}{}({}).{}'.format(plots_dir, file_name, repeated_plots, extension)
            repeated_plots += 1
    print('Saving plot in {0}'.format(save_path))
    return save_path

def plot_multiple_data_sets_overlayed(
    many_data_to_plot: list[data_id], 
    file_save_name: str,
    plotter: Callable[[plt.Axes, str, str, list], None], 
    plotter_args: list = None,
    one_time_plotter: Callable[[plt.Axes, str, str], None] = None,
    one_time_plotter_args: list = None
    ) -> None:
    plot_params.square()

    fig, ax = plt.subplots(1,1)
    for index, data_to_plot in enumerate(many_data_to_plot):
        if plotter_args == None:
            plotter(ax, data_to_plot.name, data_to_plot.legend_name)
        elif type(plotter_args[0]) == list:
            plotter(ax, data_to_plot.name, data_to_plot.legend_name, *plotter_args[index])
        else:
            plotter(ax, data_to_plot.name, data_to_plot.legend_name, *plotter_args)
    
    if one_time_plotter != None:
        if plotter_args == None:
            one_time_plotter(ax, data_to_plot.name)
        else:
            one_time_plotter(ax, data_to_plot.name, *one_time_plotter_args)
    
    plot_name = UI_helper.name_the_plot()
    fig.tight_layout()
    plt.legend(
        bbox_to_anchor = [0.45, 0.5], 
        loc='upper left',
    )

    f_name = '{}{}'.format(plot_name, file_save_name)
    save_path = define_save_plot(global_plots_dir, f_name)
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
            plotter(axs, data_to_plot.name, data_to_plot.legend_name)
        else:
            plotter(axs, data_to_plot.name, data_to_plot.legend_name, *plotter_args)
    plt.legend()
    plt.xlabel('Time [$\mathrm{T_{bin}}$]')
    
    plot_name = UI_helper.name_the_plot()

    fname = '{}disc_eccentricity_averages_{}'.format(plot_name, var)
    save_path = define_save_plot(global_plots_dir, fname)
    plt.savefig(save_path)
    plt.close()

def plot_instability_zone_for_line_plot(ax: plt.Axes, n_min: float, n_max: float) -> None:
    t_min = n_min * nts
    t_max = n_max * nts
    radius = Unit_conv.distance(instability_limit_astrophysical_object.radius)
    ax.plot([Unit_conv.time(t_min), Unit_conv.time(t_max)], [radius, radius], 'b--', linewidth = 2, label='instability limit')

def plot_Kep47b_for_line_plot(ax: plt.Axes, data_name: str, n_min, n_max) -> None:
    t_min = n_min * nts
    t_max = n_max * nts

    radius = Unit_conv.distance(Kep47b_astrophysical_object.radius)
    ax.plot([Unit_conv.time(t_min), Unit_conv.time(t_max)], [radius, radius], 'g--', linewidth = 2, label='kep47b final orbit')

def plot_instability_zone_for_twoD_sigma(ax: plt.Axes) -> None:
    radius = Unit_conv.distance(instability_limit_astrophysical_object.radius)
    theta = np.arange(0, 2*np.pi, 0.01)
    ax.plot(radius*np.cos(theta), radius*np.sin(theta), '--w', lw=1.5)

def plot_Kepler_47_planets_for_twoD_sigma(ax: plt.Axes) -> None:
    K_47b_radius = Unit_conv.distance(Kep47b_astrophysical_object.radius)
    K_47c_radius = Unit_conv.distance(Kep47c_astrophysical_object.radius)
    K_47d_radius = Unit_conv.distance(Kep47d_astrophysical_object.radius)
    theta = np.arange(0, 2*np.pi, 0.01)
    ax.plot(K_47b_radius*np.cos(theta), K_47b_radius*np.sin(theta), '#fffdcfff', lw=3)
    ax.plot(K_47d_radius*np.cos(theta), K_47d_radius*np.sin(theta), '#fffdcfff', lw=3)
    ax.plot(K_47c_radius*np.cos(theta), K_47c_radius*np.sin(theta), '#fffdcfff', lw=3)