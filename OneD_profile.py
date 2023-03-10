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

def plot_one() -> None:
    data_name = UI_helper.selectDataToPlot()
    directories = Navigation_helper.Directories(data_name)
    data_file_to_plot = UI_helper.selectDataFileToPlot(directories.out_dir)
    plot(data_name, data_file_to_plot)

def plot_many() -> None:
    data_name = UI_helper.selectDataToPlot()
    directories = Navigation_helper.Directories(data_name)
    many_data_files_to_plot = UI_helper.selectManyDataFilesToPlot(directories.out_dir)
    for data_file_to_plot in many_data_files_to_plot:
        plot(data_name, data_file_to_plot)

def animate() -> None:
    plot_params.square()
    data_name = UI_helper.selectDataToPlot()
    directories = Navigation_helper.Directories(data_name)
    n_min, n_max = UI_helper.selectAnimateRange(directories.out_dir)

    data = pluto.Pluto(directories.out_dir)

    fig = plt.figure()
    camera = Camera(fig)

    for n in range(n_min, n_max+1):  
        plot_the_data(fig, n, data_name)
        camera.snap()

    save_path = '{}1d_profile_{}_ANIMATION_{}-{}.gif'.format(directories.plots_dir, var, n_min, n_max)
    repeated_plots = 1
    while(os.path.isfile(save_path)):
        save_path = '{}1d_profile_{}_ANIMATION_{}-{}({}).gif'.format(directories.plots_dir, var, n_min, n_max, repeated_plots)
        repeated_plots += 1
    animation = camera.animate()
    animation.save(save_path)

def plot(data_name: str, data_file_to_plot: int) -> None: 
    plot_params.square()
    n = data_file_to_plot
    directories = Navigation_helper.Directories(data_name)

    fig = plt.figure()
    plot_the_data(fig, n, data_name)

    save_path = '{}1d_profile_{}_{}.png'.format(directories.plots_dir, var, n)
    repeated_plots = 1
    while(os.path.isfile(save_path)):
        save_path = '{}1d_profile_{}_{}({}).png'.format(directories.plots_dir, var, n, repeated_plots)
        repeated_plots += 1
    print('Saving plot in {0}'.format(save_path))
    plt.savefig(save_path)
    plt.close()

def plot_n_bodies(rho_max: float, data_name: str, data_file_to_plot: int) -> None:
    n = data_file_to_plot
    directories = Navigation_helper.Directories(data_name)
    num_bodies = Data_parser_helper.findNumBodies(directories.out_dir)
    for body_id in range(2, num_bodies):
        (
            time,
            a,
            e,
            omega,
            anomoly,
        ) = Data_parser_helper.getNbodyInformation_out(data_name, objects_in_kep47[body_id]) 
        plt.plot([a[n], a[n]], [0., rho_max], color='red', label=f'{objects_in_kep47[body_id].name} position')

def plot_the_data(fig: plt.Figure, n: int, data_name: str):
    directories = Navigation_helper.Directories(data_name)
    data = pluto.Pluto(directories.out_dir)
    var_data = data.primitive_variable(var, n)[0,:,:] #* data.units['density']
    print(np.argwhere(np.isnan(var_data)))
    #sigma[np.isnan(sigma)] = 1.0
    R, Phi = np.meshgrid(data.grid['faces'][0], data.grid['faces'][1])
    sig_ref = data.user_def_parameters['SIGMA_REF']
    alpha_sigma = data.user_def_parameters['ALPHA_SIGMA']

    scaled_average_azimuthal_rho = (np.mean(var_data, axis=0)/ sig_ref) / (R[0,:-1]**-alpha_sigma)
    plot_n_bodies(max(scaled_average_azimuthal_rho), data_name, n)

    if var == 'rho':
        plt.plot(R[1,:-1], scaled_average_azimuthal_rho, color='orange', label=r'scaled $\Sigma$')
        # plt.plot(R[1,:-1], np.mean(var_data, axis=0)/ sig_ref, color='orange')
        # plt.plot(R[0], (R[0]**-alpha_sigma), '-', color='blue') 
        plt.ylabel(r'$\Sigma \, [g/cm^2]$')
    elif var == 'vx1':
        plt.plot(R[1,:-1], np.mean(var_data, axis=0), color='orange')
        plt.ylabel(var)
    elif var == 'vx2':
        #TODO: figure out the units here to do this properly
        Kep_vel = 0 # tools.kepler_velocity(R[:-1,:-1][0], data) * data.read_units()['velocity']
        Kep_ang_vel = Kep_vel / R[:-1,:-1][0]
        plt.plot(R[1,:-1], np.mean(var_data, axis=0) - Kep_vel, color='orange')
        plt.ylabel(var)
    plt.grid()
    plt.xlabel(r'radius [a_bin]')
    plt.legend()
    plt.title('Kep 47')

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

