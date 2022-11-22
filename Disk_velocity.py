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

def plot_one() -> None:
    data_name = UI_helper.selectDataToPlot()
    out_dir = all_data_dir + data_name + '/out'
    data_file_to_plot = UI_helper.selectDataFileToPlot(out_dir)
    plot(data_name, data_file_to_plot)

def plot(data_name: str, data_file_to_plot: int) -> None:  
    n = data_file_to_plot
    data_parent_dir = all_data_dir + data_name
    out_dir = all_data_dir + data_name + '/out' 
    plots_dir = data_parent_dir + '/Plots/'

    data = pluto.Pluto(out_dir)
    vx1 = data.primitive_variable('vx2', n)[0,:,:] #* data.units['density']
    print(np.argwhere(np.isnan(vx1)))
    #sigma[np.isnan(sigma)] = 1.0
    R, Phi = np.meshgrid(data.grid['faces'][0], data.grid['faces'][1])
    # sig_r = data.user_def_parameters['SIGMA_REF']
    # alpha_sigma = data.user_def_parameters['ALPHA_SIGMA']

    plt.plot(R[1,:-1], np.mean(vx1, axis=0), color='orange')
    plt.xlabel(r'radius [a_bin]')
    plt.ylabel(r'$ Vx1 \, [g/cm^2]$')
    plt.title('Kep 47')

    save_path = '{}{}_vx1_profile_{}.png'.format(plots_dir, data_name, n)
    repeated_plots = 1
    while(os.path.isfile(save_path)):
        save_path = '{}{}_vx1_profile_{}({}).png'.format(plots_dir, data_name, n, repeated_plots)
        repeated_plots += 1
    plt.savefig(save_path)
    plt.close()


if __name__ == '__main__':
    plotters = [plot_one]
    func_index = UI_helper.selectFunctionsToRun(plotters)
    eval('{}()'.format(plotters[func_index].__name__))