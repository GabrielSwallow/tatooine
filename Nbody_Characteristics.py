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
import Disc_Characteristics
import plotter_helper

class possible_data_to_plot():
    eccentricity = 'eccentricity'
    semi_major_axis = 'semi major axis'


def plot_one_using_dat() -> None:
    data_name = UI_helper.selectDataToPlot()
    plot_using_dat(data_name)

def plot_one_using_dat_planet_and_cavity() -> None:
    data_name = UI_helper.selectDataToPlot()
    plot_using_dat_planet_and_cavity(data_name)

def plot_one_using_out() -> None:
    data_name = UI_helper.selectDataToPlot()
    plot_using_out(data_name)

def plot_using_dat_planet_and_cavity(data_name: str) -> None:
    directories = Navigation_helper.Directories(data_name)
    n_min, n_max = UI_helper.selectPlottingRange(directories.out_dir)
    num_avg = UI_helper.select_averaging_length()

    data_to_plot_list = [
        # possible_data_to_plot.eccentricity,
        possible_data_to_plot.semi_major_axis,
    ]

    # Jupiter_astrophysical_object.id = 2
    # Kep47b_astrophysical_object.id = 3

    objects_to_plot_list = [
        # Jupiter_astrophysical_object,
        # KepStar1_astrophysical_object,
        # KepStar2_astrophysical_object,
        Kep47b_astrophysical_object,
        # Kep47c_astrophysical_object,
        # Kep47d_astrophysical_object,
        cavity_astrophysical_object,
    ]   

    num_plots = len(data_to_plot_list)
    plot_params.one_by_N_subplots(num_plots)
    fig, axs = plt.subplots(1, num_plots, sharex= 'all')
    if num_plots == 1: axs = [axs]

    show_final_data = True
    show_instability_limit = False
    for object_to_plot in objects_to_plot_list:
        for j, data_to_plot in enumerate(data_to_plot_list):
            if object_to_plot == cavity_astrophysical_object:
                Disc_Characteristics.plot_the_data_gap_parameters_out(axs[j], data_name, 'cavity', n_min, n_max, 10, data_to_plot)
            else:
                plot_the_data_using_dat(axs[j], data_name, object_to_plot.name, object_to_plot, n_min, n_max, num_avg, data_to_plot, show_final_data, show_instability_limit)
                show_final_data = False
                show_instability_limit = False

    plt.legend()
    fig.tight_layout()
    
    save_path = '{}orbital_elements_many_objects_{}-{}.png'.format(directories.plots_dir, n_min, n_max)
    repeated_plots = 0
    while(os.path.isfile(save_path)):
        save_path = '{}orbital_elements_many_objects_{}-{}({}).png'.format(directories.plots_dir, n_min, n_max, repeated_plots)
        repeated_plots += 1
    
    print('Saving plot in {0}'.format(save_path))
    fig.savefig(save_path)
    plt.close(fig)

def plot_using_dat(data_name: str) -> None:
    directories = Navigation_helper.Directories(data_name)
    object = UI_helper.selectObjectToPlot(directories.out_dir)
    n_min, n_max = UI_helper.selectPlottingRange(directories.out_dir)
    num_avg = UI_helper.select_averaging_length()

    plot_params.one_by_two_subplot()
    fig, axs = plt.subplots(1, 2, sharex= 'all')
    plot_the_data_using_dat(axs[0], data_name, '', object, n_min, n_max, num_avg, possible_data_to_plot.eccentricity)
    plot_the_data_using_dat(axs[1], data_name, '', object, n_min, n_max, num_avg, possible_data_to_plot.semi_major_axis)
    plt.legend()
    fig.tight_layout()
    
    save_path = '{}obj{}_orbital_elements_dat_{}-{}.png'.format(directories.plots_dir, object.id, n_min, n_max)
    repeated_plots = 0
    while(os.path.isfile(save_path)):
        save_path = '{}obj{}_orbital_elements_dat_{}-{}({}).png'.format(directories.plots_dir, object.id, n_min, n_max, repeated_plots)
        repeated_plots += 1
    
    print('Saving plot in {0}'.format(save_path))
    fig.savefig(save_path)
    plt.close(fig)

def plot_the_data_using_dat(
        ax: plt.Axes, 
        data_name: str, 
        legend_name: str,
        object: astrophysical_object,
        n_min: int, 
        n_max: int,
        num_avg: int,
        data_to_plot: str,
        show_final_data: bool = False,
        show_instability_limit: bool = False,
    ):

    t_min = n_min * nts
    t_max = n_max * nts
    (
        time,
        a,
        e,
        period,
        mass,
    ) = Data_parser_helper.getNbodyInformation_dat(data_name, object.id) 

    i_min, i_max = tools.time_split(time, t_min, t_max)

    if data_to_plot == 'eccentricity':
        e_split_rolling_average, t_split_rolling_average = tools.rolling_average(num_avg, e[i_min:i_max], time)
        ax.plot(t_split_rolling_average, e_split_rolling_average, label = legend_name)
        ax.set_ylabel(r'Eccentricity')
    elif data_to_plot == 'semi major axis':
        a_split_rolling_average, t_split_rolling_average = tools.rolling_average(num_avg, a[i_min:i_max], time)
        ax.plot(t_split_rolling_average, a_split_rolling_average, label = legend_name)
        if show_final_data:
            plotter_helper.plot_Kep47b_for_line_plot(ax, t_split_rolling_average[0], t_split_rolling_average[-1])
        if show_instability_limit:
            plotter_helper.plot_instability_zone_for_line_plot(ax, t_split_rolling_average[0], t_split_rolling_average[-1])
        ax.set_ylabel(r'Semi-Major Axis [$\mathrm{a_{bin}}$]')
    ax.set_xlabel('Time [$\mathrm{T_{bin}}$]')

    return object, n_min, n_max

def plot_using_out(data_name: str) -> None:
    plot_params.square()
    directories = Navigation_helper.Directories(data_name)

    object = UI_helper.selectObjectToPlot(directories.out_dir)
    n_min, n_max = UI_helper.selectPlottingRange(directories.out_dir)
    
    (
        time,
        a,
        e,
        omega,
        anomaly,
    ) = Data_parser_helper.getNbodyInformation_out(directories.data_name, object) 
    
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
    fig.suptitle('Orbital Elements of {}'.format(object.name))
    
    fig.tight_layout()
    
    save_path = '{}obj{}_orbital_elements_{}-{}.png'.format(directories.plots_dir, object.id, n_min, n_max)
    repeated_plots = 0
    while(os.path.isfile(save_path)):
        save_path = '{}obj{}_orbital_elements_{}-{}({}).png'.format(directories.plots_dir, n_min, n_max, repeated_plots)
        repeated_plots += 1
    
    print('Saving plot in {0}'.format(save_path))
    fig.savefig(save_path)
    plt.close(fig)

def plot_resonance_dat() -> None:
    plot_params.square()
    data_name = UI_helper.selectDataToPlot()
    directories = Navigation_helper.Directories(data_name)

    obj_list, obj_des_list = UI_helper.selectObjectsToPlot(directories.out_dir)
    n_min, n_max = UI_helper.selectPlottingRange(directories.out_dir)
    t_min = n_min * nts
    t_max = n_max * nts

    (
        time_0,
        a_0,
        e_0,
        period_0,
        mass_0,
    ) = Data_parser_helper.getNbodyInformation_dat(data_name, obj_list[0])
    
    (
        time_1,
        a_1,
        e_1,
        period_1,
        mass_1,
    ) = Data_parser_helper.getNbodyInformation_dat(data_name, obj_list[1])

    a_norm = a_1/a_0
    period_norm = a_norm ** 3/2
    e_norm = e_1/e_0

    i_min, i_max = tools.time_split(time_0, t_min, t_max)

    fig, axs = plt.subplots(2,1, sharex = 'all')

    axs[0].plot(time_0[i_min:i_max], period_norm[i_min:i_max])
    axs[0].set_title(r'Period ratio')

    axs[1].plot(time_0[i_min:i_max], a_norm[i_min:i_max])
    axs[1].set_title(r'Semi-Major Axis ratio')

    axs[0].grid()
    axs[1].grid()

    #axs[1].plot(time_0[i_min:i_max], e_norm[i_min:i_max])
    #axs[1].set_title(r'Eccentricity ratio [$\mathrm{e}$]')
    axs[1].set(xlabel = r'Time [$\mathrm{T_{bin}}$]')
    fig.suptitle('Resonances of kepler-47{}'.format(obj_des_list[1]))

    fig.tight_layout()

    save_path = '{}obj{}_obj{}_resonances_dat_{}-{}.png'.format(directories.plots_dir, obj_list[0], obj_list[1], n_min, n_max)
    repeated_plots = 0
    while(os.path.isfile(save_path)):
        save_path = '{}obj{}_obj{}_resonances_dat_{}-{}({}).png'.format(directories.plots_dir, obj_list[0], obj_list[1], n_min, n_max, repeated_plots)
        repeated_plots += 1
    
    print('Saving plot in {0}'.format(save_path))
    fig.savefig(save_path)
    plt.close(fig)

def plot_resonance_dat_fit() -> None:
    plot_params.square()
    data_name = UI_helper.selectDataToPlot()
    directories = Navigation_helper.Directories(data_name)

    obj_list, obj_des_list = UI_helper.selectObjectsToPlot(directories.out_dir)
    n_min, n_max = UI_helper.selectPlottingRange(directories.out_dir)
    t_min = n_min * nts
    t_max = n_max * nts

    (
        time_0,
        a_0,
        e_0,
        period_0,
        mass_0,
    ) = Data_parser_helper.getNbodyInformation_dat(data_name, obj_list[0])
    
    (
        time_1,
        a_1,
        e_1,
        period_1,
        mass_1,
    ) = Data_parser_helper.getNbodyInformation_dat(data_name, obj_list[1])

    a_norm = a_1/a_0
    avg = UI_helper.select_averaging_length()
    a_norm = tools.roll_avg(a_norm, avg)
    time_norm = tools.roll_avg(time_0, avg)
    period_norm = a_norm ** 3/2
    e_norm = e_1/e_0

    a_fit = np.polyfit(time_norm, a_norm, 1)
    period_fit = np.polyfit(time_norm, period_norm, 1)


    i_min, i_max = tools.time_split(time_norm, t_min, t_max)

    fig, axs = plt.subplots(2,1, sharex = 'all')

    axs[0].plot(time_norm[i_min:i_max], period_norm[i_min:i_max])
    #axs[0].plot(time_norm[i_min:i_max], period_fit[0] * time_norm[i_min:i_max] + period_fit[1], 'r--')
    #period_fit_text = axs[0].text(0.95, 0.95, r'fit={}'.format(period_fit[1]),
    #        color='k',
    #        size='x-large',
    #        ha='right',
    #        va='top',transform=axs[0].transAxes)
    axs[0].set_title(r'Period ratio')

    axs[1].plot(time_norm[i_min:i_max], a_norm[i_min:i_max])
    #axs[1].plot(time_norm[i_min:i_max], a_fit[0] * time_norm[i_min:i_max] + a_fit[1], 'r--')
    #a_fit_text = axs[1].text(0.95, 0.95, r'fit={}'.format(a_fit[1]),
    #        color='k',
    #        size='x-large',
    #        ha='right',
    #        va='top',transform=axs[1].transAxes)
    axs[1].set_title(r'Semi-Major Axis ratio')

    axs[0].grid()
    axs[1].grid()

    #axs[1].plot(time_0[i_min:i_max], e_norm[i_min:i_max])
    #axs[1].set_title(r'Eccentricity ratio [$\mathrm{e}$]')
    axs[1].set(xlabel = r'Time [$\mathrm{T_{bin}}$]')
    fig.suptitle('Resonances of kepler-47{}'.format(obj_des_list[1]))

    fig.tight_layout()

    save_path = '{}obj{}_obj{}_resonances_dat_{}-{}.png'.format(directories.plots_dir, obj_list[0], obj_list[1], n_min, n_max)
    repeated_plots = 0
    while(os.path.isfile(save_path)):
        save_path = '{}obj{}_obj{}_resonances_dat_{}-{}({}).png'.format(directories.plots_dir, obj_list[0], obj_list[1], n_min, n_max, repeated_plots)
        repeated_plots += 1
    
    print('Saving plot in {0}'.format(save_path))
    fig.savefig(save_path)
    plt.close(fig)
    
if __name__ == '__main__':
    plotters = [
        plot_one_using_dat, 
        plot_one_using_dat_planet_and_cavity,
        plot_one_using_out, 
        plot_resonance_dat, 
        plot_resonance_dat_fit]
    func_index = UI_helper.selectFunctionsToRun(plotters)
    eval('{}()'.format(plotters[func_index].__name__))