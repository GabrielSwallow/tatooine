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

def calculate_Torque(data, data_index: int, data_name: str, obj: int = 2):
    n = data_index
    directories = Navigation_helper.Directories(data_name)

    sigma = data.primitive_variable('rho', n)[0,:,:] #* data.units['density']
    sig_r = data.user_def_parameters['SIGMA_REF']
    alpha_sigma = data.user_def_parameters['ALPHA_SIGMA']

    R, Phi = np.meshgrid(data.grid['faces'][0], data.grid['faces'][1])
    x = tools.x_coord(R, Phi)
    y = tools.y_coord(R, Phi)
    dV = data.grid['dV']

    (
        _,
        a_nbody_list,
        _,
        _,
        anomoly_nbody_list,
    ) = Data_parser_helper.getNbodyInformation_out(directories.out_dir, obj)

    a_nbody = a_nbody_list[n]
    anomoly_nbody = anomoly_nbody_list[n]
    x_nbody = tools.x_coord(a_nbody, anomoly_nbody)
    y_nbody = tools.y_coord(a_nbody, anomoly_nbody)

    nbody_mass = Data_parser_helper.get_planet_mass(data_name)[obj-2]

    norm = (sigma*dV*nbody_mass) / ((R - a_nbody)**2)
    t1 = y * (x - x_nbody)
    t2 = x * (y - y_nbody)
    torque_grid = norm * (t1 + t2)