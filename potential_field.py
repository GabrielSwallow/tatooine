import pluto
import numpy as np
import copy
from tools import x_coord, y_coord, constants
from numba import jit 
import UI_helper
import Navigation_helper

import pickle
import os

G = constants['G']

def get_x_and_y_grids(x0_len, x1_len, r_grid, phi_grid):
    x_grid = np.zeros([x0_len, x1_len], float)
    y_grid = np.zeros([x0_len, x1_len], float)
    for i in range(x0_len):
        for j in range(x1_len):
            x_grid[i][j] = x_coord(r_grid[i][j], phi_grid[i][j])
            y_grid[i][j] = y_coord(r_grid[i][j], phi_grid[i][j])
    return x_grid, y_grid

@jit
def calc_potential_due_to_cell(
        x_pos_test_particle: float, 
        y_pos_test_particle: float, 
        x_pos_cell: float,
        y_pos_cell: float,
        mass_cell: float,
        ):
    return G * mass_cell / (( (x_pos_test_particle-x_pos_cell)**2 + (y_pos_test_particle-y_pos_cell)**2 )**0.5)

@jit
def calc_potential_at_pos(i, j, x0_len, x1_len, x_pos_test_particle, y_pos_test_particle, x_grid, y_grid, mass_grid):
    potential = 0
    for m in range(x0_len):
        for n in range(x1_len):
            if m==i or n==j:
                continue
            x_pos_cell = x_grid[m][n]
            y_pos_cell = y_grid[m][n]
            mass_cell = mass_grid[m][n]
            calc_potential_due_to_cell(x_pos_test_particle, y_pos_test_particle, x_pos_cell,y_pos_cell, mass_cell)
            # potential += G * mass_grid[m][n] / (( (x_pos_test_particle-x_pos_cell)**2 + (y_pos_test_particle-y_pos_cell)**2 )**0.5)
    return potential

@jit
def calc_row_of_potentials(i, x0_len, x1_len, x_grid, y_grid, mass_grid):
    potential_field_row = np.zeros(x1_len)
    for j in range(x1_len):
            x_pos_test_particle = x_grid[i][j]
            y_pos_test_particle = y_grid[i][j]
            potential_field_row[j] = calc_potential_at_pos(i, j, x0_len, x1_len, x_pos_test_particle, y_pos_test_particle, x_grid, y_grid, mass_grid)
    return potential_field_row

def calculate_potential_field(): 
    data_name = UI_helper.selectDataToPlot()
    directories = Navigation_helper.Directories(data_name)
    out_file_number = UI_helper.selectDataFileToPlot(directories.out_dir)
    data = pluto.Pluto(directories.out_dir)
    pickle_name = directories.extra_data_dir + 'potential_field_{}.pkl'.format(out_file_number)

    if os.path.isfile(pickle_name):
        print('\npickle already exists\n')
        return 

    var_data = data.primitive_variable('rho', out_file_number)[0,:,:] #* data.units['density']
    x0_len = len(var_data)
    x1_len = len(var_data[0])

    r_grid, phi_grid = np.meshgrid(data.grid['centers'].X1, data.grid['centers'].X2)
    x_grid, y_grid = get_x_and_y_grids(x0_len, x1_len, r_grid, phi_grid)
    potential_field = np.zeros([x0_len, x1_len], float)

    dVx1, dVx2 = np.meshgrid(data.grid['dV'][0], data.grid['dV'][1])
    dV = dVx1*dVx2

    mass_grid = var_data*dV

    for i in range(x0_len):
        print(i)
        potential_field[i] = calc_row_of_potentials(i, x0_len, x1_len, x_grid, y_grid, mass_grid)

    with open(pickle_name, 'wb') as f:
        pickle.dump(potential_field, f)
    return potential_field

if __name__ == '__main__':
    x = calculate_potential_field()