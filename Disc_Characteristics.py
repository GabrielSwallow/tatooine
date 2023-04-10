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
import Gap_Finder as gap

def plot_disc_e_avg() -> None:
    plot_params.square()
    data_name = UI_helper.selectDataToPlot()
    directories = Navigation_helper.Directories(data_name)

    fig, axs = plt.subplots(1,1)
    plot_the_data_disc_e_avg(axs, data_name)
    fig.tight_layout()

    f_name = 'disc_eccentricity_averages'
    save_path = plotter_helper.define_save_plot(directories.plots_dir, f_name)
    fig.savefig(save_path)
    plt.close(fig)

def plot_one_multiple_data_sets_overlayed_disc_e_avg() -> None:
    plot_params.square()
    data_ids = UI_helper.select_many_data_ids_to_overlay(choose_out_file = False)
    avg_num = UI_helper.select_averaging_length()
    plotter_args = [avg_num]
    plotter_helper.plot_multiple_data_sets_overlayed(data_ids, 'disc_e_avg', plot_the_data_disc_e_avg, plotter_args)

def plot_the_data_disc_e_avg(ax: plt.Axes, data_name: str, legend_name: str = '', avg_num: int = 1) -> None:
    directories = Navigation_helper.Directories(data_name)
    n_min, n_max = UI_helper.selectPlottingRange(directories.out_dir)
    t_min = n_min * nts
    t_max = n_max * nts

    (time, _, e, _, _, _, _, _, _, _) = Data_parser_helper.get_averages_data(data_name)
    i_min, i_max = tools.time_split(time, t_min, t_max)
    e_split_rolling_average, t_split_rolling_average = tools.rolling_average(avg_num, e[i_min:i_max], time[i_min:i_max])
    ax.plot(Unit_conv.time(t_split_rolling_average), e_split_rolling_average, label = f'{legend_name}')
    ax.set_ylabel('Eccentricity')
    ax.set_xlabel('Time [' + Unit_conv.time_label() + ']')

def plot_gap_parameters_out() -> None:
    plot_params.two_by_one_subplot()
    data_name = UI_helper.selectDataToPlot()
    directories = Navigation_helper.Directories(data_name)
    fig, axs = plt.subplots(2,1, sharex = 'all')

    n_min, n_max = UI_helper.selectPlottingRange(directories.out_dir)
    avg_num = UI_helper.select_averaging_length()
    plot_the_data_gap_parameters_out(axs[0], data_name, '', n_min, n_max, avg_num, data_to_plot='semi major axis')
    plot_the_data_gap_parameters_out(axs[1], data_name, '', n_min, n_max, avg_num, data_to_plot='eccentricity')
    axs[1].set(xlabel = 'Time [$\mathrm{T_{bin}}$]')

    fig.tight_layout()
    # fig.suptitle('Disc Parameters of Kepler-47')

    fname = 'gap_eccentricity_{}-{}'.format(n_min, n_max)
    save_path = plotter_helper.define_save_plot(directories.plots_dir, fname)
    fig.savefig(save_path)
    plt.close(fig)

def plot_one_multiple_data_sets_gap_parameters_out() -> None:
    data_ids = UI_helper.select_many_data_ids_to_overlay(choose_out_file = False)
    n_min, n_max = UI_helper.selectPlottingRange()
    avg_num = UI_helper.select_averaging_length()

    plotter_args = [n_min, n_max, avg_num]
    plotter_helper.plot_multiple_data_sets_overlayed(data_ids, 'gap_parameters_out', plot_the_data_gap_parameters_out, plotter_args)

def plot_the_data_gap_parameters_out(
        ax: plt.Axes, 
        data_name: str, 
        legend_name: str = '',
        n_min: int = 0,
        n_max: int = 0,
        avg_num: int = 1,
        data_to_plot: str = 'eccentricity'
    ) -> tuple[int]:
    directories = Navigation_helper.Directories(data_name)
    
    e = []
    a = []
    t = []
    data = pluto.Pluto(directories.out_dir)
    for n in range(n_min, n_max+1):
        t.append(n*nts)
        var_data = data.primitive_variable(var, n)[0,:,:]
        R, Phi = np.meshgrid(data.grid['faces'][0], data.grid['faces'][1])
        (a_n, e_n, _, _, _, _, _) = gap.ellipse_find(R, Phi, var_data)
        e.append(abs(e_n))
        a.append(abs(a_n))
    pass
    if data_to_plot == 'eccentricity':
        rolling_avg_e, rolling_avg_time = tools.rolling_average(avg_num, np.array(e), np.array(t))
        ax.plot(Unit_conv.time(rolling_avg_time), rolling_avg_e, label=legend_name)
        ax.set_ylabel('Eccentricity')
    elif data_to_plot == 'semi major axis':
        rolling_avg_a, rolling_avg_time = tools.rolling_average(avg_num, a, t)
        ax.plot(Unit_conv.time(rolling_avg_time), Unit_conv.distance(rolling_avg_a), label=legend_name)
        ax.set_ylabel('Semi-Major Axis [' + Unit_conv.distance_label() + ']')
    ax.set_xlabel('Time [' + Unit_conv.time_label() + ']')
    return n_min, n_max

def calculate_total_disc_accretion() -> float:
    data_name = UI_helper.selectDataToPlot()
    directories = Navigation_helper.Directories(data_name)
    (time, _, _, _, _, _, _, _, _, accs) = Data_parser_helper.get_averages_data(data_name)
    total_disc_accretion = [accs[0][i] + accs[1][i] + accs[2][i] for i in range(len(accs[0]))]
    acc_tot = np.trapz(total_disc_accretion, time)
    # print('total accretion onto disc = {}'.format(acc_tot))
    print_statement = 'average accretion rate for {} = {} {} / {}'.format(
        data_name,
        Unit_conv.mass(acc_tot, 'Earth') / Unit_conv.time(time[-1], 'years'), 
        Unit_conv.mass_label('Earth'),
        Unit_conv.time_label('years')
        )
    print(print_statement)
    return acc_tot, print_statement

def plot_difference_between_min_and_max_rho_same_radius():
    plot_params.square()
    data_name = UI_helper.selectDataToPlot()
    directories = Navigation_helper.Directories(data_name)
    out_file_number = UI_helper.selectDataFileToPlot(directories.out_dir)

    data = pluto.Pluto(directories.out_dir)
    rho = data.primitive_variable("rho", out_file_number)[0,:,:]
    num_radii_points = len(rho[0])
    r_vals = data.grid['centers'].X1

    phi_min = []
    phi_max = []
    # phi_diff = []
    
    for r_index in range(num_radii_points):
        data_at_fixed_r = rho[:, r_index]
        phi_min_val = min(data_at_fixed_r)
        phi_max_val = max(data_at_fixed_r)
        phi_min.append(phi_min_val)
        phi_max.append(phi_max_val)
        # phi_diff.append(phi_max_val - phi_min_val)

    fig = plt.figure()    

    fig.tight_layout()

    plt.plot(
        tools.Unit_conv.distance(r_vals), 
        tools.Unit_conv.surface_density(np.array(phi_min), conv_unit_mass='grams', conv_unit_distance='cm'), 
        label='min surface density'
    )
    plt.plot(
        tools.Unit_conv.distance(r_vals), 
        tools.Unit_conv.surface_density(np.array(phi_max), conv_unit_mass='grams', conv_unit_distance='cm'), 
        label='max surface density'
    )
    # plt.plot(r_vals, phi_diff)
    plt.legend()
    plt.yscale('log')
    plt.xlabel('radius [{}]'.format(tools.Unit_conv.distance_label()))
    plt.ylabel('surface density [{}]'.format(tools.Unit_conv.surface_density_label(conv_unit_mass='grams', conv_unit_distance='cm')))


    fname = 'disc_eccentricity_vs_radius_{}'.format(out_file_number)
    save_path = plotter_helper.define_save_plot(directories.plots_dir, fname, 'png')
    fig.savefig(save_path)
    plt.close(fig)


if __name__ == '__main__':
    plotters = [
        plot_disc_e_avg, 
        plot_one_multiple_data_sets_overlayed_disc_e_avg, 
        plot_gap_parameters_out, 
        plot_one_multiple_data_sets_gap_parameters_out,
        calculate_total_disc_accretion,
        plot_difference_between_min_and_max_rho_same_radius,
    ]
    func_index = UI_helper.selectFunctionsToRun(plotters)
    eval('{}()'.format(plotters[func_index].__name__))