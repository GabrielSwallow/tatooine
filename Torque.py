import pluto
import tools
import numpy as np
import os 
import matplotlib.pyplot as plt

import UI_helper
import Navigation_helper
import Data_parser_helper
from Global_variables import *
import plot_params

def print_torque_from_calculation():
    data_name = UI_helper.selectDataToPlot()
    directories = Navigation_helper.Directories(data_name)
    data_file_to_plot = UI_helper.selectDataFileToPlot(directories.out_dir)
    obj_index, _ = UI_helper.selectObjectToPlot(directories.out_dir)
    return calculate_torque(data_name, data_file_to_plot, obj_index)

def plot_torque_from_averages():
    plot_params.square()
    data_name = UI_helper.selectDataToPlot()
    directories = Navigation_helper.Directories(data_name)
    obj_index, _ = UI_helper.selectObjectToPlot(directories.out_dir)
    avg_num = UI_helper.select_averaging_length()
    (time, _, _,  _, _, _, _, _, torque_list, _) = Data_parser_helper.get_averages_data(data_name)
    obj_torque = torque_list[obj_index-2]
    rolling_average_torque = tools.rolling_average(obj_torque, avg_num)

    fig = plt.figure()
    plt.plot(time, abs(rolling_average_torque), label='torque')
    plt.xlabel('Time (binary orbits)')
    plt.ylabel('Torque (unknown units)')
    plt.yscale('log')
    plt.legend()
    plt.grid()

    save_path = '{}{}_obj{}_torque.png'.format(directories.plots_dir, data_name, obj_index)
    repeated_plots = 0
    while(os.path.isfile(save_path)):
        save_path = '{}{}_obj{}_torque({}).png'.format(directories.plots_dir, data_name, obj_index, repeated_plots)
        repeated_plots += 1
    
    print('Saving plot in {0}'.format(save_path))
    fig.savefig(save_path)
    plt.close(fig)

def plot_torque_from_averages_with_inner_and_outer():
    plot_params.square()
    data_name = UI_helper.selectDataToPlot()
    directories = Navigation_helper.Directories(data_name)
    obj_index, _ = UI_helper.selectObjectToPlot(directories.out_dir)
    (time, _, _,  _, _, _, _, _, torque_list, _) = Data_parser_helper.get_averages_data(data_name)
    obj_torque_inner = torque_list[(obj_index-2)*2]
    obj_torque_outer = torque_list[(obj_index-2)*2 + 1]

    fig = plt.figure()
    plt.plot(time, abs(obj_torque_inner), label='torque')
    plt.plot(time, abs(obj_torque_outer), label='torque')
    plt.xlabel('Time (binary orbits)')
    plt.ylabel('Torque (unknown units)')
    plt.yscale('log')
    plt.legend()
    plt.grid()

    save_path = '{}{}_obj{}_torque.png'.format(directories.plots_dir, data_name, obj_index)
    repeated_plots = 0
    while(os.path.isfile(save_path)):
        save_path = '{}{}_obj{}_torque({}).png'.format(directories.plots_dir, data_name, obj_index, repeated_plots)
        repeated_plots += 1
    
    print('Saving plot in {0}'.format(save_path))
    fig.savefig(save_path)
    plt.close(fig)

def plot_torque_from_calculation():
    plot_params.square()
    data_name = UI_helper.selectDataToPlot()
    directories = Navigation_helper.Directories(data_name)
    obj_index, _ = UI_helper.selectObjectToPlot(directories.out_dir)

    inter_torque_list = []
    outer_torque_list = []
    time_list = []
    i = 47
    f = 50
    for n in range(Navigation_helper.findMaxFileNumber(directories.out_dir)):
        inner_torque, outer_torque, time = calculate_torque(data_name, n, obj_index)
        inter_torque_list.append(abs(inner_torque))
        outer_torque_list.append(abs(outer_torque))
        time_list.append(time)
    fig = plt.figure()
    
    plt.plot(time_list, inter_torque_list, label='inner torque')
    plt.plot(time_list, outer_torque_list, label='outer torque')
    plt.xlabel('Time (binary orbits)')
    plt.ylabel('Torque (unknown units)')
    plt.yscale('log')
    plt.legend()
    plt.grid()

    save_path = '{}{}_obj{}_torque.png'.format(directories.plots_dir, data_name, obj_index)
    repeated_plots = 0
    while(os.path.isfile(save_path)):
        save_path = '{}{}_obj{}_torque({}).png'.format(directories.plots_dir, data_name, obj_index, repeated_plots)
        repeated_plots += 1
    
    print('Saving plot in {0}'.format(save_path))
    fig.savefig(save_path)
    plt.close(fig)

def calculate_torque(data_name: str, data_index: int, obj_index: int = 2):
    directories = Navigation_helper.Directories(data_name)
    data = pluto.Pluto(directories.out_dir)
    n = data_index

    sigma = data.primitive_variable('rho', n)[0,:,:] #* data.units['density']
    # sig_r = data.user_def_parameters['SIGMA_REF']
    # alpha_sigma = data.user_def_parameters['ALPHA_SIGMA']

    R, Phi = np.meshgrid(data.grid['faces'][0], data.grid['faces'][1])
    volume_grid_R, volume_grid_phi = np.meshgrid(data.grid['dV'][0], data.grid['dV'][1])
    x = tools.x_coord(R, Phi)[:-1, :-1]
    y = tools.y_coord(R, Phi)[:-1, :-1]

    (
        _,
        x_nbody_list,
        y_nbody_list,
        z,
        vx,
        vy,
        vz,
    ) = Data_parser_helper.getNbodyCoordinates(data_name, obj_index)

    (
        dbl_num_orbits_per_out, 
        analysis_num_orbits_per_log,
    ) = Data_parser_helper.getAnalysisOutMetaInfo(data_name)
    analysis_logs_per_dbl_out = int(dbl_num_orbits_per_out / analysis_num_orbits_per_log)
    time = dbl_num_orbits_per_out * n

    x_nbody = x_nbody_list[analysis_logs_per_dbl_out * n]
    y_nbody = y_nbody_list[analysis_logs_per_dbl_out * n]
    r_nbody = tools.r_coord(x_nbody, y_nbody)

    (_, _, _, _, nbody_mass_list)  = Data_parser_helper.getNbodyInformation_dat(data_name, obj_index)
    if not isinstance(nbody_mass_list, np.ndarray):
        print("\nWARNING: No mass column in nbody_orbital_elements.dat: using planet's initial mass\n")
        nbody_mass = Data_parser_helper.get_initial_planet_masses(data_name)[obj_index-2]
    else:
        nbody_mass = nbody_mass_list[analysis_logs_per_dbl_out * n]


    norm = (sigma*volume_grid_R*nbody_mass) / ((R[:-1, :-1] - r_nbody)**2)
    t1 = y * (x - x_nbody)
    t2 = x * (y - y_nbody)
    torque_grid = norm * (t1 + t2)

    outer_torque = 0.0
    inner_torque = 0.0
    for i, torque_row in enumerate(torque_grid):
        for j, torque_cell in enumerate(torque_row):
            if R[i][j] < r_nbody:
                inner_torque += torque_cell
            else:
                outer_torque += torque_cell
    return inner_torque, outer_torque, time

def get_average_disk_velocity_at_r(data_name: str, data_index: int, wavenumber: int) -> list[float]:
    directorites = Navigation_helper.Directories(data_name)
    data = pluto.Pluto(directorites.out_dir)
    var_data = data.primitive_variable('vx2', data_index)[0,:,:] #* data.units['density']
    return [np.mean(var_data[:, r]) for r in range(684)]

def get_planet_velocity(data_name: str, data_index: int, wavenumber: int, obj_index: int = 2) -> float:
    n = data_index
    (
        _,
        x_nbody_list,
        y_nbody_list,
        z,
        vx,
        vy,
        vz,
    ) = Data_parser_helper.getNbodyCoordinates(data_name, obj_index)

    (
        dbl_num_orbits_per_out, 
        analysis_num_orbits_per_log,
    ) = Data_parser_helper.getAnalysisOutMetaInfo(data_name)
    analysis_logs_per_dbl_out = int(dbl_num_orbits_per_out / analysis_num_orbits_per_log)
    time = dbl_num_orbits_per_out * n

    x_nbody = x_nbody_list[analysis_logs_per_dbl_out * n]
    y_nbody = y_nbody_list[analysis_logs_per_dbl_out * n]
    r_nbody = tools.r_coord(x_nbody, y_nbody)
    return tools.kepler_velocity(r_nbody)

def calculate_Lindblad_torque(data_name: str, data_index: int, wavenumber: int, obj_index: int = 2):
    disk_velocities = get_average_disk_velocity_at_r(data_name, data_index, wavenumber)
    planet_velocity = get_planet_velocity(data_name, data_index, wavenumber, obj_index)
    print(planet_velocity)
    plt.plot(range(len(disk_velocities)), disk_velocities)
    plt.show()


if __name__ == '__main__':
    plotters = [
        plot_torque_from_averages, 
        plot_torque_from_averages_with_inner_and_outer, 
        plot_torque_from_calculation
    ]
    func_index = UI_helper.selectFunctionsToRun(plotters)
    eval('{}()'.format(plotters[func_index].__name__))

    # inner_torque, outer_torque, time = get_torque_from_out_file()
    # print(inner_torque, outer_torque, time)
    

    # calculate_Lindblad_torque(
    #     data_name='accretion_with_in_large',
    #     data_index=50,
    #     wavenumber=2,
    #     obj_index=2
    # )

    # plot_torque()