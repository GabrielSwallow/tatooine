import pluto
import tools
import numpy as np
import os 

import UI_helper
import Navigation_helper
import Data_parser_helper
from Global_variables import *

def torque_from_data():
    data_name = UI_helper.selectDataToPlot()
    directories = Navigation_helper.Directories(data_name)
    data_file_to_plot = UI_helper.selectDataFileToPlot(directories.out_dir)
    obj_index, _ = UI_helper.selectObjectToPlot(directories.out_dir)
    return calculate_Torque(data_name, data_file_to_plot, obj_index)

def calculate_Torque(data_name: str, data_index: int, obj_index: int = 2):
    directories = Navigation_helper.Directories(data_name)
    data = pluto.Pluto(directories.out_dir)
    n = data_index

    sigma = data.primitive_variable('rho', n)[0,:,:] #* data.units['density']
    sig_r = data.user_def_parameters['SIGMA_REF']
    alpha_sigma = data.user_def_parameters['ALPHA_SIGMA']

    R, Phi = np.meshgrid(data.grid['faces'][0], data.grid['faces'][1])
    volume_grid_R, volume_grid_phi = np.meshgrid(data.grid['dV'][0], data.grid['dV'][1])
    x = tools.x_coord(R, Phi)[:-1, :-1]
    y = tools.y_coord(R, Phi)[:-1, :-1]
    dV = data.grid['dV']

    (
        _,
        a_nbody_list,
        _,
        _,
        anomoly_nbody_list,
    ) = Data_parser_helper.getNbodyInformation_out(data_name, obj_index)
    # TODO: use nbody_coordinates.dat instead for exact pos

    a_nbody = a_nbody_list[n]
    anomoly_nbody = anomoly_nbody_list[n]
    x_nbody = tools.x_coord(a_nbody, anomoly_nbody)
    y_nbody = tools.y_coord(a_nbody, anomoly_nbody)

    nbody_mass = Data_parser_helper.get_planet_masses(data_name)[obj_index-2]

    norm = (sigma*volume_grid_R*nbody_mass) / ((R[:-1, :-1] - a_nbody)**2)
    t1 = y * (x - x_nbody)
    t2 = x * (y - y_nbody)
    torque_grid = norm * (t1 + t2)

    outer_torque = 0.
    inner_torque = 0.
    for i, torque_row in enumerate(torque_grid):
        for j, torque_cell in enumerate(torque_row):
            if R[i][j] < a_nbody:
                inner_torque += torque_cell
            else:
                outer_torque += torque_cell
    return inner_torque, outer_torque

if __name__ == '__main__':
    inner_torque, outer_torque = torque_from_data()
    print(inner_torque, outer_torque)