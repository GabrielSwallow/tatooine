#!/usr/bin/env python3

import pluto
import tools
import numpy as np
import matplotlib.pyplot as plt
import UI_helper
import Navigation_helper
import os 
from celluloid import Camera

all_data_dir = '../Data/'

def plot_one():
    data_name = UI_helper.selectDataToPlot()
    out_dir = all_data_dir + data_name + '/out'
    data_file_to_plot = UI_helper.selectDataFileToPlot(out_dir)
    plot(data_name, data_file_to_plot)

def plot_many():
    data_name = UI_helper.selectDataToPlot()
    out_dir = all_data_dir + data_name + '/out'
    many_data_files_to_plot = UI_helper.selectManyDataFilesToPlot(out_dir)
    for d in many_data_files_to_plot:
        plot(data_name, many_data_files_to_plot[d])

def animate():
    data_name = UI_helper.selectDataToPlot()
    out_dir = all_data_dir + data_name + '/out'
    n_max = Navigation_helper.findMaxFileName(out_dir)

    data_parent_dir = all_data_dir + data_name
    out_dir = all_data_dir + data_name + '/out' 
    plots_dir = data_parent_dir + '/Plots/'
    data = pluto.Pluto(out_dir)

    fig = plt.figure()
    camera = Camera(fig)

    for n in range(n_max):  
        
        sigma = data.primitive_variable('rho', n)[0,:,:] #* data.units['density']
        print(np.argwhere(np.isnan(sigma)))
        #sigma[np.isnan(sigma)] = 1.0
        R, Phi = np.meshgrid(data.grid['faces'][0], data.grid['faces'][1])
        sig_r = data.user_def_parameters['SIGMA_REF']

        plt.plot(R[1,:-1], np.mean(sigma, axis=0)/ sig_r, color='orange')
        plt.plot(R[0], R[0]**-0.5, '-', color='blue')
        #plt.ylim(0,0.1)
        print(np.max(np.mean(sigma, axis=0)/ data.units['density'] ))
        plt.xlim(0,20)
        plt.xlabel(r'radius [a_bin]')
        plt.ylabel(r'$\Sigma \, [g/cm^2]$')
        plt.title('Kep 47')
        camera.snap()

    save_path = '{}{}_1d_profile_ANIMATION.gif'.format(plots_dir, data_name)
    repeated_plots = 1
    while(os.path.isfile(save_path)):
        save_path = '{}{}_1d_profile_ANIMATION({}).gif'.format(plots_dir, data_name, repeated_plots)
        repeated_plots += 1
    animation = camera.animate()
    animation.save(save_path)

def plot(data_name: str, data_file_to_plot: int) -> None:  
    n = data_file_to_plot
    data_parent_dir = all_data_dir + data_name
    out_dir = all_data_dir + data_name + '/out' 
    plots_dir = data_parent_dir + '/Plots/'   

    data = pluto.Pluto(out_dir)
    sigma = data.primitive_variable('rho', n)[0,:,:] #* data.units['density']
    print(np.argwhere(np.isnan(sigma)))
    #sigma[np.isnan(sigma)] = 1.0
    R, Phi = np.meshgrid(data.grid['faces'][0], data.grid['faces'][1])
    sig_r = data.user_def_parameters['SIGMA_REF']

    fig = plt.figure()
    plt.figure()
    plt.plot(R[1,:-1], np.mean(sigma, axis=0)/ sig_r  )
    plt.plot(R, R**-0.5, '-')
    #plt.ylim(0,0.1)
    print(np.max(np.mean(sigma, axis=0)/ data.units['density'] ))
    plt.xlim(0,20)
    plt.xlabel(r'radius [a_bin]')
    plt.ylabel(r'$\Sigma \, [g/cm^2]$')
    plt.title('Kep 47')

    save_path = '{}{}_1d_profile_{}.png'.format(plots_dir, data_name, n)
    repeated_plots = 1
    while(os.path.isfile(save_path)):
        save_path = '{}{}_1d_profile_{}({}).png'.format(plots_dir, data_name, n, repeated_plots)
        repeated_plots += 1

    fig.savefig(save_path)


if __name__ == '__main__':
    plotters = [plot_one, plot_many, animate]
    func_index = UI_helper.selectFunctionsToRun(plotters)
    eval('{}()'.format(plotters[func_index].__name__))





# plt.figure()
# v_phi = data.primitive_variable('vx2', n)[0,:,:]
# plt.plot(R[1,:-1], np.mean(v_phi, axis=0))
# #plt.ylim(0,1)
# plt.xlim(0,20)
# plt.xlabel(r'radius [a_bin]')
# plt.ylabel(r'$v_{\phi} [T_{bin}]$')
# plt.title('')
# plt.show()

#%%

# for n in range(0,15):
#     print(n)
#     sigma = data.primitive_variable('rho', n)[0,:,:] #* data.units['density']
#     print(np.argwhere(np.isnan(sigma)))
#     #sigma[np.isnan(sigma)] = 1.0
#     R, Phi = np.meshgrid(data.grid['faces'][0], data.grid['faces'][1])
#     sig_r = data.user_def_parameters['SIGMA_REF']

#     plt.figure()
#     plt.plot(R[1,:-1], np.mean(sigma, axis=0)/ sig_r  )
#     plt.plot(R, R**-0.5, '-')
#     #plt.ylim(0,0.1)
#     print(np.max(np.mean(sigma, axis=0)/ data.units['density'] ))
#     plt.xlim(0,20)
#     plt.xlabel(r'radius [a_bin]')
#     plt.ylabel(r'$\Sigma \, [g/cm^2]$')
#     plt.title('Kep 47')
#     plt.show()

