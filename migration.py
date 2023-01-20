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
    delta_t_list = np.array([time[i+1] - time[i] for i in range(len(a)-1)])
    a_dot_list = delta_a_list / delta_t_list
    a_dot_over_a_list = a_dot_list/a[:-1]
    plt.plot(time[::10][:-1], a_dot_list[::10])
    plt.show()

if __name__ == '__main__':
    calculate_migration()