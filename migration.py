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

def calculate_migration():
    data_name = UI_helper.selectDataToPlot()
    directories = Navigation_helper.Directories(data_name)
    obj_index, _ = UI_helper.selectObjectToPlot(directories.out_dir)


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

    new_list_length = int(len(delta_a_list)/10)
    time_list = np.array([time[i*10] for i in range(new_list_length)])
    larger_delta_a_list = np.array([sum(delta_a_list[i*10:(i*10)+10]) for i in range(new_list_length)])
    larger_a_delta_t_list = np.array([sum(a_delta_t_list[i*10:(i*10)+10]) for i in range(new_list_length)])

    a_dot_over_a_list = larger_delta_a_list / larger_a_delta_t_list

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
    calculate_migration()