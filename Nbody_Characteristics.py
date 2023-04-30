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
    objects_to_plot_list = UI_helper.select_object_config_to_plot() 
    data_to_plot_list = UI_helper.select_a_or_e_to_plot()
    plot_using_dat_planet_and_cavity(data_name, objects_to_plot_list, data_to_plot_list)

def plot_one_using_out() -> None:
    data_name = UI_helper.selectDataToPlot()
    plot_using_out(data_name)

def plot_many_data_id_using_dat() -> None:
    data_ids = UI_helper.select_many_data_ids_to_overlay()
    first_data_point_directories = Navigation_helper.Directories(data_ids[0].name)
    object = UI_helper.selectObjectToPlot(first_data_point_directories.out_dir)
    num_avg = UI_helper.select_averaging_length()

    show_instability_limit = UI_helper.show_instability_limit()
    show_47b_final_orbit = UI_helper.show_47b_final_orbit()

    data_to_plot_list = UI_helper.select_a_or_e_to_plot()
    if len(data_to_plot_list) == 2: raise Exception('plot_many_data_id_using_dat only takes a or e, not both at the moment')
    if data_to_plot_list[0] == possible_data_to_plot.eccentricity: data_name_short = 'e'
    elif data_to_plot_list[0] == possible_data_to_plot.semi_major_axis: data_name_short = 'a'

    plotter_args = []
    for data in data_ids:
        directories = Navigation_helper.Directories(data.name)
        n_min, n_max = 0, Navigation_helper.findMaxFileNumber(directories.out_dir)
        plotter_args.append([object, n_min, n_max, num_avg, data_to_plot_list[0], show_47b_final_orbit, show_instability_limit, data.plot_colour])

    
    Kep47b_for_line_plot_args = [n_min, n_max]
    
    if data_to_plot_list[0] == possible_data_to_plot.eccentricity or not show_47b_final_orbit:
        plotter_helper.plot_multiple_data_sets_overlayed(
            data_ids, 
            '{}_{}_evolution'.format(object.name, data_name_short), 
            plot_the_data_using_dat, 
            plotter_args, 
        )
    else:
        plotter_helper.plot_multiple_data_sets_overlayed(
            data_ids, 
            '{}_{}_evolution'.format(object.name, data_name_short), 
            plot_the_data_using_dat, 
            plotter_args, 
            plotter_helper.plot_Kep47b_for_line_plot,
            Kep47b_for_line_plot_args,
        )

def plot_using_dat_planet_and_cavity(data_name: str, objects_to_plot_list: list[astrophysical_object], data_to_plot_list: list[str]) -> None:
    directories = Navigation_helper.Directories(data_name)
    n_min, n_max = UI_helper.selectPlottingRange(directories.out_dir)
    num_avg = UI_helper.select_averaging_length()
    show_instability_limit = UI_helper.show_instability_limit()
    show_47b_final_orbit = UI_helper.show_47b_final_orbit()

    num_plots = len(data_to_plot_list)
    plot_params.one_by_N_subplots(num_plots)
    fig, axs = plt.subplots(1, num_plots, sharex= 'all')
    if num_plots == 1: axs = [axs]

    for object_to_plot in objects_to_plot_list:
        for j, data_to_plot in enumerate(data_to_plot_list):
            if object_to_plot == cavity_astrophysical_object:
                Disc_Characteristics.plot_the_data_gap_parameters_out(axs[j], data_name, 'cavity', n_min, n_max, 10, data_to_plot)
            else:
                plot_the_data_using_dat(axs[j], data_name, object_to_plot.name, object_to_plot, n_min, n_max, num_avg, data_to_plot, show_47b_final_orbit, show_instability_limit)
                show_47b_final_orbit = False
                show_instability_limit = False

    plt.legend()
    fig.tight_layout()

    if len(data_to_plot_list) == 2: suffix = 'a_and_e'
    if data_to_plot_list[0] == possible_data_to_plot.eccentricity: suffix = 'e'
    elif data_to_plot_list[0] == possible_data_to_plot.semi_major_axis: suffix = 'a'

    
    fname = '{}_many_objects_{}-{}'.format(suffix, n_min, n_max)
    save_path = plotter_helper.define_save_plot(directories.plots_dir, fname)
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
    

    fname = 'obj{}_orbital_elements_dat_{}-{}'.format(object.id, n_min, n_max)
    save_path = plotter_helper.define_save_plot(directories.plots_dir, fname)
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
        e_split_rolling_average, t_split_rolling_average = tools.rolling_average(num_avg, e[i_min:i_max], time[i_min:i_max])
        ax.plot(Unit_conv.time(t_split_rolling_average), e_split_rolling_average, label = legend_name)
        ax.set_ylabel(r'Eccentricity')
    elif data_to_plot == 'semi major axis':
        a_split_rolling_average, t_split_rolling_average = tools.rolling_average(num_avg, a[i_min:i_max], time[i_min:i_max])
        ax.plot(Unit_conv.time(t_split_rolling_average), Unit_conv.distance(a_split_rolling_average), label = legend_name)
        if show_final_data:
            plotter_helper.plot_Kep47b_for_line_plot(ax, data_name, n_min, n_max)
        if show_instability_limit:
            plotter_helper.plot_instability_zone_for_line_plot(ax, t_split_rolling_average[0], t_split_rolling_average[-1])
        ax.set_ylabel('Semi-Major Axis [' + Unit_conv.distance_label() + ']')
        # ax.set_ylim([0,15])
    ax.set_xlabel('Time [' + Unit_conv.time_label() + ']')

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
    
    fname = 'obj{}_orbital_elements_{}-{}'.format(object.id, n_min, n_max)
    save_path = plotter_helper.define_save_plot(directories.plots_dir, fname)
    fig.savefig(save_path)
    plt.close(fig)

def plot_resonance_dat() -> None:
    plot_params.square()
    data_name = UI_helper.selectDataToPlot()
    directories = Navigation_helper.Directories(data_name)
    obj_list, obj_des_list = UI_helper.selectObjectsToPlot(directories.out_dir)
    n_min, n_max = UI_helper.selectPlottingRange(directories.out_dir)
    avg_num = UI_helper.select_averaging_length()

    fig, ax = plt.subplots(1,1, sharex = 'all')
    plot_the_data_resonance_dat(ax, data_name, '', obj_list[0], obj_list[1], n_min, n_max, avg_num)
    fig.tight_layout()
    fname = 'obj{}_obj{}_resonances_dat_{}-{}'.format(obj_list[0], obj_list[1], n_min, n_max)
    save_path = plotter_helper.define_save_plot(directories.plots_dir, fname)
    fig.savefig(save_path)
    plt.close(fig) 

def plot_the_data_resonance_dat(
    ax: plt.Axes, 
    data_name: str, 
    legend_name: str = '', 
    obj_id_1: int = 2, 
    obj_id_2: int = 3, 
    n_min: int = 0, 
    n_max: int = 0,
    avg_num: int = 1,
    ) -> None:
    directories = Navigation_helper.Directories(data_name)
    t_min = n_min * nts
    t_max = n_max * nts

    (
        time_0,
        a_0,
        _,
        _,
        _,
    ) = Data_parser_helper.getNbodyInformation_dat(data_name, obj_id_1)
    
    (
        time_1,
        a_1,
        _,
        _,
        _,
    ) = Data_parser_helper.getNbodyInformation_dat(data_name, obj_id_2)

    a_norm = a_1/a_0
    i_min, i_max = tools.time_split(time_0, t_min, t_max)
    a_norm, time_norm = tools.rolling_average(avg_num, a_norm[i_min:i_max], time_0[i_min:i_max])
    period_norm = np.power(a_norm, 3/2)

    period_norm = [p if p>1 else 1/p for p in period_norm]

    ax.plot(Unit_conv.time(time_norm), period_norm)
    ax.set_ylabel('ratio of periods')
    ax.set_xlabel('Time [' + Unit_conv.time_label() + ']')
    # ax.set_xlim([200, 240])
    # ax.set_ylim([0,3])


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

    fname = 'obj{}_obj{}_resonances_dat_{}-{}'.format(obj_list[0], obj_list[1], n_min, n_max)
    save_path = plotter_helper.define_save_plot(directories.plots_dir, fname)
    fig.savefig(save_path)
    plt.close(fig)
    
if __name__ == '__main__':
    plotters = [
        plot_one_using_dat, 
        plot_one_using_dat_planet_and_cavity,
        plot_one_using_out, 
        plot_resonance_dat, 
        plot_resonance_dat_fit,
        plot_many_data_id_using_dat]
    func_index = UI_helper.selectFunctionsToRun(plotters)
    eval('{}()'.format(plotters[func_index].__name__))