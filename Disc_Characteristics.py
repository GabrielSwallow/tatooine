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
    data_name = UI_helper.selectDataToPlot()
    directories = Navigation_helper.Directories(data_name)

    fig = plt.figure()
    plot_the_data_disc_e_avg(fig, data_name)
    fig.tight_layout()

    save_path = '{}disc_eccentricity_averages.png'.format(directories.plots_dir)
    repeated_plots = 0
    while (os.path.isfile(save_path)):
        save_path = '{}disc_eccentricity_averages({}).png'.format(directories.plots_dir, repeated_plots)
        repeated_plots += 1
    
    print('Saving plot in {0}'.format(save_path))
    fig.savefig(save_path)
    plt.close(fig)

def plot_one_multiple_data_sets_overlayed_disc_e_avg() -> None:
    data_ids = UI_helper.select_many_data_ids_to_overlay(choose_out_file = False)
    plotter_helper.plot_multiple_data_sets_overlayed(data_ids, plot_the_data_disc_e_avg)

# def plot_multiple_data_sets_overlayed_disc_e_avg(many_data_to_plot: list[data_id]) -> None:
#     plot_params.square()

#     fig = plt.figure()
#     for data_to_plot in many_data_to_plot:
#         plot_the_data_disc_e_avg(fig, data_to_plot.name, data_to_plot.legend_name)
    
#     plot_name = UI_helper.name_the_plot() + '_'

#     save_path = '{}disc_eccentricity_averages_{}.png'.format(global_plots_dir+plot_name, var)
#     repeated_plots = 1
#     while(os.path.isfile(save_path)):
#         save_path = '{}disc_eccentricity_averages_{}({}).png'.format(global_plots_dir+plot_name, var, repeated_plots)
#         repeated_plots += 1
#     print('Saving plot in {0}'.format(save_path))
#     plt.legend()
#     plt.ylabel('gap eccentricity')
#     plt.xlabel('Time [$\mathrm{T_{bin}}$]')
#     plt.savefig(save_path)
#     plt.close()

def plot_the_data_disc_e_avg(fig: plt.Figure, data_name: str, legend_name: str = '') -> None:
    plot_params.square()
    directories = Navigation_helper.Directories(data_name)
    n_min, n_max = UI_helper.selectPlottingRange(directories.out_dir)
    t_min = n_min * nts
    t_max = n_max * nts

    (time, _, e, _, _, _, _, _, _, _) = Data_parser_helper.get_averages_data(data_name)
    i_min, i_max = tools.time_split(time, t_min, t_max)
    plt.plot(time[i_min:i_max], e[i_min:i_max], label = f'{legend_name} gap eccentricity')


def plot_gap_parameters_out() -> None:
    data_name = UI_helper.selectDataToPlot()
    directories = Navigation_helper.Directories(data_name)
    fig, axs = plt.subplots(2,1, sharex = 'all')
    n_min, n_max = plot_the_data_gap_parameters_out(fig, axs, data_name)

    fig.tight_layout()

    save_path = '{}gap_eccentricity_{}-{}.png'.format(directories.plots_dir, n_min, n_max)
    repeated_plots = 0
    while (os.path.isfile(save_path)):
        save_path = '{}gap_eccentricity_{}-{}({}).png'.format(directories.plots_dir, n_min, n_max, repeated_plots)
        repeated_plots += 1
    
    print('Saving plot in {0}'.format(save_path))
    fig.savefig(save_path)
    plt.close(fig)

def plot_one_multiple_data_sets_gap_parameters_out() -> None:
    data_ids = UI_helper.select_many_data_ids_to_overlay(choose_out_file = False)
    plotter_helper.plot_multiple_data_sets_overlayed_subplots(data_ids, plot_the_data_gap_parameters_out)

def plot_the_data_gap_parameters_out(fig: plt.Figure, axs: plt.Axes, data_name: str, legend_name: str = '') -> tuple[int]:
    directories = Navigation_helper.Directories(data_name)
    n_min, n_max = UI_helper.selectPlottingRange(directories.out_dir)
    avg_num = UI_helper.select_averaging_length()

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
    
    axs[0].plot(t[:-avg_num], tools.rolling_average(a, avg_num), label=legend_name)
    axs[0].set_title('Semi-Major Axis [$\mathrm{a_{bin}}$]')

    axs[1].plot(t[:-avg_num], tools.rolling_average(e, avg_num), label=legend_name)
    axs[1].set_title('Eccentricity [$\mathrm{e}$]')

    axs[1].set(xlabel = 'Time [$\mathrm{T_{bin}}$]')
    fig.suptitle('Disc Parameters of Kepler-47')

    return n_min, n_max

if __name__ == '__main__':
    plotters = [
        plot_disc_e_avg, 
        plot_one_multiple_data_sets_overlayed_disc_e_avg, 
        plot_gap_parameters_out, 
        plot_one_multiple_data_sets_gap_parameters_out
    ]
    func_index = UI_helper.selectFunctionsToRun(plotters)
    eval('{}()'.format(plotters[func_index].__name__))