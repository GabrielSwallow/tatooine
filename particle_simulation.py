import numpy as np
from scipy import spatial
from scipy.integrate import odeint

import UI_helper
import Navigation_helper
import pluto
from potential_field import load_potential_field
from numba import jit

def combine_two_matricies(M1, M2):
    x0_len = len(M1)
    x1_len = len(M1[0])
    return np.array([np.array([ np.array([M1[i][j], M2[i][j]]) for j in range(x1_len)]) for i in range(x0_len)])
     

def calculate_force_field(data: pluto.Pluto, potential_field:list[list[float]]) -> list[list[list[float]]]:
    r_vals = data.grid['centers'].X1
    phi_vals = data.grid['centers'].X2
    r_grid, phi_grid = np.meshgrid(r_vals, phi_vals)

    x0_len = len(r_grid)
    x1_len = len(r_grid[0])

    dVdr = np.array([np.gradient(potential_field[i], r_grid[i]) for i in range(x0_len)])
    dVdphi_transposed = np.array([np.gradient(potential_field[:, i], phi_grid[:, i]) for i in range(x1_len)])
    dVdphi = np.matrix.transpose(dVdphi_transposed)
    dVdphi_grad_term = np.array([
        np.array([dVdphi[i][j] * (1/r_grid[i][j])
        for j in range(x1_len)]) 
        for i in range(x0_len)])
    
    return combine_two_matricies(dVdr, dVdphi_grad_term)

def find_nearest_cell(data: pluto.Pluto, particle_position: list[float]):
    r_vals = data.grid['centers'].X1
    phi_vals = data.grid['centers'].X2
    r_grid, phi_grid = np.meshgrid(r_vals, phi_vals)
    grid_of_coords = combine_two_matricies(r_grid, phi_grid)
    flattened_grid = np.matrix.flatten(grid_of_coords)
    flattened_coords = np.array( [[flattened_grid[2*i], flattened_grid[2*i+1]] for i in range(int(len(flattened_grid)/2))] )

    x0_len = len(grid_of_coords)
    x1_len = len(grid_of_coords[0])

    tree = spatial.KDTree(flattened_coords)
    query = tree.query(particle_position)
    print('distance to point = ', query[0])
    numerical_coord = flattened_coords[query[1]]
    index_of_coord_i = int(query[1]/x0_len)
    index_of_coord_j = query[1]%x1_len
    return numerical_coord, [index_of_coord_i, index_of_coord_j]

def calculate_acceleration(data: pluto.Pluto, force_field, particle_position, particle_mass):
    numerical_coord, index_of_coord = (data, particle_position)
    force = force_field[index_of_coord[0], index_of_coord[1]]
    return force/particle_mass


def propagate_particle(data: pluto.Pluto, force_field, particle_mass):
    def ode_to_solve(vars, t):
        vx,vy,x,y = vars
        if y<0:
            y += 2*np.pi
        elif y>2*np.pi:
            y -= 2*np.pi
        return *calculate_acceleration(data, force_field, [x,y], particle_mass), vx, vy
    # ode_to_solve = lambda vx,vy,x,y: calculate_acceleration(data, force_field, [x,y], particle_mass)
    # position = lambda x,y:
    solution = odeint(
        func = ode_to_solve,
        y0 = [10.,0.,0.,3.],
        t=np.linspace(0,100,1000)
    )
    return solution


if __name__ == '__main__':
    data_name = UI_helper.selectDataToPlot()
    directories = Navigation_helper.Directories(data_name)
    data_file_name = UI_helper.select_potential_field_to_plot(directories.data_name)
    potential_field = load_potential_field(data_name, data_file_name)
    data = pluto.Pluto(directories.out_dir)

    particle_position = [1,1]
    force_field = calculate_force_field(data, potential_field)
    find_nearest_cell(data, particle_position)
    sol = propagate_particle(data, force_field, particle_mass=1)