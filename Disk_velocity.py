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

def plot_one() -> None:
    data_name = UI_helper.selectDataToPlot()
    out_dir = all_data_dir + data_name + '/out'
    data_file_to_plot = UI_helper.selectDataFileToPlot(out_dir)
    plot(data_name, data_file_to_plot)

def plot(data_name: str, data_file_to_plot: int) -> None:  
    plot_params.square()
    n = data_file_to_plot
    directories = Navigation_helper.Directories(data_name)
    # data_parent_dir = all_data_dir + data_name
    # out_dir = all_data_dir + data_name + '/out' 
    # plots_dir = data_parent_dir + '/Plots/'

    data = pluto.Pluto(directories.out_dir)
    vx1 = data.primitive_variable('vx2', n)[0,:,:] #* data.units['density']
    print(np.argwhere(np.isnan(vx1)))
    #sigma[np.isnan(sigma)] = 1.0
    R, Phi = np.meshgrid(data.grid['faces'][0], data.grid['faces'][1])
    # sig_r = data.user_def_parameters['SIGMA_REF']
    # alpha_sigma = data.user_def_parameters['ALPHA_SIGMA']
    a_max = 40
    R_cells = 684
    a_to_plot = 15
    R_cell_to_plot = int(a_to_plot * (R_cells/a_max))
    # plt.plot(R[1,:-1], np.mean(vx1, axis=0), color='orange')
    # plt.plot(R[1,:-1], vx1[0], color='orange')
    plt.plot(Phi[:-1, R_cell_to_plot], vx1[:, R_cell_to_plot], color='orange')
    # plt.plot(R[1,:-1], np.mean(vx1, axis=0), color='orange')
    plt.xlabel(r'radius [a_bin]')
    plt.ylabel(r'$ Vx1 \, [g/cm^2]$')
    plt.title('Kep 47')

    fname = '{}vx1_profile_{}'.format(directories.plots_dir, n)
    save_path = plotter_helper.define_save_plot(fname)
    plt.savefig(save_path)
    plt.close()

def animate() -> None:
    plot_params.square()
    data_name = UI_helper.selectDataToPlot()
    directories = Navigation_helper.Directories(data_name)
    n_max = Navigation_helper.findMaxFileNumber(directories.out_dir)

    max_file = Navigation_helper.findMaxFileNumber(directories.out_dir)

    fig = plt.figure()
    camera = Camera(fig)

    for n in range(max_file+1):
        data = pluto.Pluto(directories.out_dir)
        vx1 = data.primitive_variable('vx2', n)[0,:,:] #* data.units['density']
        print(np.argwhere(np.isnan(vx1)))
        #sigma[np.isnan(sigma)] = 1.0
        R, Phi = np.meshgrid(data.grid['faces'][0], data.grid['faces'][1])
        # sig_r = data.user_def_parameters['SIGMA_REF']
        # alpha_sigma = data.user_def_parameters['ALPHA_SIGMA']

        a_max = 40
        R_cells = 684
        a_to_plot = 15
        R_cell_to_plot = int(a_to_plot * (R_cells/a_max))
        # plt.plot(R[1,:-1], np.mean(vx1, axis=0), color='orange')
        # plt.plot(R[1,:-1], vx1[0], color='orange')
        plt.plot(Phi[:-1, R_cell_to_plot], vx1[:, R_cell_to_plot], color='orange')
        plt.xlabel(r'radius [a_bin]')
        plt.ylabel(r'$ Vx1 \, [g/cm^2]$')
        plt.title('Kep 47')
        camera.snap()

    fname = '{}vx1_profile_{}'.format(directories.plots_dir, n)
    save_path = plotter_helper.define_save_plot(fname, 'gif')
    animation = camera.animate()
    animation.save(save_path)

if __name__ == '__main__':
    plotters = [plot_one, animate]
    func_index = UI_helper.selectFunctionsToRun(plotters)
    eval('{}()'.format(plotters[func_index].__name__))