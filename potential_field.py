import pluto
import numpy as np
import copy
from tools import x_coord, y_coord, constants, binary
from numba import jit 
import UI_helper
import Navigation_helper
import Data_parser_helper

import seaborn
import matplotlib.pyplot as plt
import matplotlib.colors as colors
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
def calc_potential_due_to_binary(
    x_pos_test_particle: float, 
    y_pos_test_particle: float, 
    ):
    return G * 1 / (( (x_pos_test_particle)**2 + (y_pos_test_particle)**2 )**0.5)

@jit
def calc_potential_due_to_body(
    x_pos_test_particle: float, 
    y_pos_test_particle: float, 
    x_pos_planet: float,
    y_pos_planet: float,
    planet_mass: float,
    ) -> float:
    return G * planet_mass / (( (x_pos_test_particle-x_pos_planet)**2 + (y_pos_test_particle-y_pos_planet)**2 )**0.5)


@jit
def calc_potential_due_to_disc_at_pos(i, j, x0_len, x1_len, x_pos_test_particle, y_pos_test_particle, x_grid, y_grid, mass_grid):
    potential = 0
    for m in range(x0_len):
        for n in range(x1_len):
            if m==i or n==j:
                continue
            x_pos_cell = x_grid[m][n]
            y_pos_cell = y_grid[m][n]
            mass_cell = mass_grid[m][n]
            potential += calc_potential_due_to_cell(x_pos_test_particle, y_pos_test_particle, x_pos_cell,y_pos_cell, mass_cell)
    # potential += calc_potential_due_to_binary(x_pos_test_particle, y_pos_test_particle)
    return potential

@jit
def calc_row_of_potentials_due_to_disc(i, x0_len, x1_len, x_grid, y_grid, mass_grid):
    potential_field_row = np.zeros(x1_len)
    for j in range(x1_len):
        x_pos_test_particle = x_grid[i][j]
        y_pos_test_particle = y_grid[i][j]
        potential_field_row[j] = calc_potential_due_to_disc_at_pos(i, j, x0_len, x1_len, x_pos_test_particle, y_pos_test_particle, x_grid, y_grid, mass_grid)
    return potential_field_row

def calculate_potential_field_due_to_disc(): 
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
        potential_field[i] = calc_row_of_potentials_due_to_disc(i, x0_len, x1_len, x_grid, y_grid, mass_grid)

    with open(pickle_name, 'wb') as f:
        print('saving data to ', pickle_name)
        pickle.dump(potential_field, f)
    return potential_field

def add_potential_due_to_nbodies():
    data_name = UI_helper.selectDataToPlot()
    directories = Navigation_helper.Directories(data_name)
    original_potential_data_file_name = UI_helper.select_potential_field_to_plot(directories.data_name)
    out_file_number = UI_helper.selectDataFileToPlot(directories.out_dir)
    add_potential_due_to_nbodies_given_info(data_name, original_potential_data_file_name, out_file_number)

def add_potential_due_to_nbodies_over_many_t():
    data_name = UI_helper.selectDataToPlot()
    directories = Navigation_helper.Directories(data_name)
    original_potential_data_file_name = UI_helper.select_potential_field_to_plot(directories.data_name)
    out_file_numbers = UI_helper.selectManyDataFilesToPlot(directories.out_dir)

    for out_file_number in out_file_numbers:
        add_potential_due_to_nbodies_given_info(data_name, original_potential_data_file_name, out_file_number)


def add_potential_due_to_nbodies_given_info(data_name: str, original_potential_data_file_name: str, out_file_number: str):
    directories = Navigation_helper.Directories(data_name)
    data = pluto.Pluto(directories.out_dir)
    disc_pickle_name = directories.extra_data_dir + original_potential_data_file_name #directories.extra_data_dir + 'potential_field_{}.pkl'.format(out_file_number)
    
    num_bodies = Data_parser_helper.findNumBodies(directories.out_dir)
    x_bodies = []
    y_bodies = []
    for i in range(num_bodies):
        _, x, y, _, _, _, _ = Data_parser_helper.getNbodyCoordinates(data_name, i)
        x_bodies.append(x)
        y_bodies.append(y)
    mass_bodies = [binary['m1'], binary['m2']]
    mass_planets = Data_parser_helper.get_initial_planet_masses(data_name)
    mass_bodies.append(*mass_planets)

    if not os.path.isfile(disc_pickle_name):
        print('\nno pickle for disc potential\n')
        return 
    
    var_data = data.primitive_variable('rho', out_file_number)[0,:,:] #* data.units['density']
    x0_len = len(var_data)
    x1_len = len(var_data[0])

    r_grid, phi_grid = np.meshgrid(data.grid['centers'].X1, data.grid['centers'].X2)
    x_grid, y_grid = get_x_and_y_grids(x0_len, x1_len, r_grid, phi_grid)

    dbl_num_orbits_per_out, analysis_num_orbits_per_log = Data_parser_helper.getAnalysisOutMetaInfo(data_name)
    num_orbits = out_file_number * dbl_num_orbits_per_out
    analysis_index = int(num_orbits / analysis_num_orbits_per_log)

    with open(directories.extra_data_dir + original_potential_data_file_name, 'rb') as f:
        potential_field = pickle.load(f)
        for i in range(x0_len):
            for j in range(x1_len):
                x_pos_test_particle = x_grid[i][j]
                y_pos_test_particle = y_grid[i][j]
                for n in range(num_bodies):
                    potential_field[i][j]+= calc_potential_due_to_body(
                        x_pos_test_particle, 
                        y_pos_test_particle, 
                        x_bodies[n][analysis_index], 
                        y_bodies[n][analysis_index], 
                        mass_bodies[n]
                    )
    
    add_binary_pickle_name = directories.extra_data_dir + 'potential_field_bodies_{}.pkl'.format(out_file_number)

    with open(add_binary_pickle_name, 'wb') as f:
        print('saving data to ', add_binary_pickle_name)
        pickle.dump(potential_field, f)
    return potential_field  

def load_potential_field(data_name: str, data_file_name: str):
    directories = Navigation_helper.Directories(data_name)
    with open(directories.extra_data_dir + data_file_name, 'rb') as f:
        potential_field = pickle.load(f)
        return potential_field

def plot_potential_field():
    data_name = UI_helper.selectDataToPlot()
    directories = Navigation_helper.Directories(data_name)
    data_file_name = UI_helper.select_potential_field_to_plot(directories.data_name)
    data = pluto.Pluto(directories.out_dir)

    potential_field = load_potential_field(data_name, data_file_name)
    var_data = data.primitive_variable('rho', 0)[0,:,:] #* data.units['density']
    r_vals = data.grid['centers'].X1
    phi_vals = data.grid['centers'].X2
    r_grid, phi_grid = np.meshgrid(r_vals, phi_vals)

    fig = plt.figure()
    # potential_field_plot = plt.imshow(phi_grid, r_grid, potential_field, origin='lower', label='potential field')
    potential_field_plot = plt.pcolormesh(phi_grid, r_grid, potential_field)
    contours_plot = plt.contour(phi_grid, r_grid, potential_field, colors='red')
    colour_bar = plt.colorbar(potential_field_plot)
    colour_bar.add_lines(contours_plot)
    # plt.subplot(projection="polar")
    # plt.ylim(1,4)

    save_path = '{}PLOT_{}.png'.format(directories.plots_dir, data_file_name)
    repeated_plots = 1
    while(os.path.isfile(save_path)):
        save_path = '{}PLOT_{}({}).png'.format(directories.plots_dir, data_file_name, repeated_plots)
        repeated_plots += 1
    print('Saving plot in {0}'.format(save_path))
    fig.savefig(save_path)
    plt.close(fig)

# def find_line_of_constant_potential(potential_field, potential):
#     radii_indices = []
#     radii_indices.append(
#         np.argmin([abs(potential_field[0][j]-potential) for j in range(len(potential_field[0]))])
#     )
#     max_delta_r = 10
#     for i in range(1, len(potential_field)):
#         min_r_range = max(0, radii_indices[-1]-max_delta_r)
#         max_r_range = min(684, radii_indices[-1]+max_delta_r)
#         radii_indices.append(
#             min_r_range+np.argmin([abs(potential_field[i][j]-potential) for j in range(min_r_range, max_r_range)])
#             )
#     return radii_indices

if __name__ == '__main__':
    functions = [
        calculate_potential_field_due_to_disc, 
        add_potential_due_to_nbodies, 
        add_potential_due_to_nbodies_over_many_t, 
        plot_potential_field
    ]
    func_index = UI_helper.selectFunctionsToRun(functions)
    eval('{}()'.format(functions[func_index].__name__))