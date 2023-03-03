import pluto
import tools
import numpy as np
import os 
import matplotlib.pyplot as plt

import UI_helper
import Navigation_helper
import Data_parser_helper
from Global_variables import *
import Torque
import plot_params

def plot_migration_one():
    data_name = UI_helper.selectDataToPlot()
    directories = Navigation_helper.Directories(data_name)
    obj_index, _ = UI_helper.selectObjectToPlot(directories.out_dir)
    avg_num = UI_helper.select_averaging_length()
    calculate_migration(data_name, obj_index, avg_num)

def plot_migration_many_avgs():
    data_name = UI_helper.selectDataToPlot()
    directories = Navigation_helper.Directories(data_name)
    obj_index, _ = UI_helper.selectObjectToPlot(directories.out_dir)
    for i in [i*5 for i in range(1, 10)]:
        calculate_migration(data_name, obj_index, i)

def calculate_migration(data_name: str, obj_index: int, avg_num: int):
    plot_params.square()

    directories = Navigation_helper.Directories(data_name)
    (
        time,
        a,
        e,
        anomoly,
        mass,
    ) = Data_parser_helper.getNbodyInformation_dat(data_name, obj_index)

    delta_a_list = np.array([a[i+1] - a[i] for i in range(len(a)-1)])
    a_delta_t_list = np.array([a[i]*(time[i+1] - time[i]) for i in range(len(a)-1)])

    list_len = int(len(delta_a_list))
    larger_delta_a_list = tools.rolling_average(delta_a_list, avg_num)
    larger_a_delta_t_list = tools.rolling_average(a_delta_t_list, avg_num)

    a_dot_over_a_list = larger_delta_a_list / larger_a_delta_t_list
    time_list = np.array(time[0:(list_len-avg_num)])

    fig = plt.figure()
    plt.plot(time_list, a_dot_over_a_list)
    plt.xlabel('Time (binary orbits)')
    plt.ylabel('a_dot/a (1/binary orbit time)')
    plt.legend()
    plt.grid()
    fig.tight_layout()
    
    save_path = '{}{}_obj{}_migration.png'.format(directories.plots_dir, data_name, obj_index)
    repeated_plots = 0
    while(os.path.isfile(save_path)):
        save_path = '{}{}_obj{}_migration({}).png'.format(directories.plots_dir, data_name, obj_index, repeated_plots)
        repeated_plots += 1
    
    print('Saving plot in {0}'.format(save_path))
    fig.savefig(save_path)
    plt.close(fig)

if __name__ == '__main__':
    plotters = [plot_migration_one, plot_migration_many_avgs]
    func_index = UI_helper.selectFunctionsToRun(plotters)
    eval('{}()'.format(plotters[func_index].__name__))