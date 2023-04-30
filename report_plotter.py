### report_plotter ###

import pluto
import tools
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib as mpl
import matplotlib.colors as colors
import argparse as argp
import os
import UI_helper
import Navigation_helper
from Global_variables import *
import plot_params as plotp
import TwoD_sigma
import OneD_profile
import Nbody_Characteristics
import Accretion
import migration
import Torque
import Data_parser_helper
import Gap_Finder as gap
import mass as m

a_crit = 2.31
kepb = 3.53
yrs = m.t['yr']['kep47']

def small_planet_typeI_blocked():

    name_data = UI_helper.selectDataToPlot()
    dir = all_data_dir + name_data + '/out'
    data = pluto.Pluto(dir)
    planet = UI_helper.selectObjectToPlot(dir)
    (
        time_dat,
        a_dat,
        e_dat,
        period_dat,
        mass_dat
    ) = Data_parser_helper.getNbodyInformation_dat(name_data, planet.id)
    _ = UI_helper.selectPlottingRange(dir)
    n_min, n_max = _[0],_[1]
    nts = int(input('nts \n'))
    t_min, t_max = nts * n_min, nts * n_max
    i_min, i_max = tools.time_split(time_dat, t_min, t_max)

    (
        gap_time,
        inner_gap_e,
        inner_gap_a,
        peak_gap_e,
        peak_gap_a,
    ) = (
        [],
        [],
        [],
        [],
        [],
    )
    for i in range(n_min, n_max+1):
        gap_time.append(i*nts*yrs)
        var_data = data.primitive_variable(var, i)[0,:,:]
        R, Phi = np.meshgrid(data.grid['faces'][0], data.grid['faces'][1])
        (
            inner_a_j,
            inner_e_j,
            _, _, _, _, _,
        ) = gap.ellipse_find(R, Phi, var_data)
        inner_gap_e.append(abs(inner_e_j))
        inner_gap_a.append(abs(inner_a_j))
        (
            peak_a_j,
            peak_e_j,
            _, _, _, _, _,
        ) = gap.ellipse_find_alt(R, Phi, var_data)
        peak_gap_e.append(abs(peak_e_j))
        peak_gap_a.append(abs(peak_a_j))

    time_dat *= yrs
    t_min *= yrs
    t_max *= yrs

    fig, ax = plt.subplots(1,1)

    l1 = ax.plot(time_dat[i_min:i_max],a_dat[i_min:i_max], 'k-')
    l2 = ax.plot(gap_time, inner_gap_a, 'b:')
    l3 = ax.plot(gap_time, peak_gap_a, 'b--')
    l4 = ax.plot([t_min, t_max], [kepb, kepb], 'g-')
    l5 = ax.plot([t_min, t_max], [a_crit, a_crit], 'r-')
    ax.set_xlabel(r'Time (years)', size =14)
    ax.set_ylabel(r'Semi-Major Axis ($a_{KEP47}$)', size = 14)
    ax.legend([r'Small Planet', r'Inner Cavity', r'Density Peak', r'Kepler 47 b', r'Critical Radius'])
    ax.grid()
    save_path = '../Report_Plot/SmallPlanetFailure_a.png'
    repeated_plots = 0
    while (os.path.isfile(save_path)):
        save_path = '../Report_Plot/SmallPlanetFailure_a({}).png'.format(repeated_plots)
        repeated_plots += 1
    print('saving plot in {}'.format(save_path))
    fig.savefig(save_path, bbox_inches = 'tight', dpi = 320)
    plt.close(fig)

    fig, ax = plt.subplots(1,1)

    l1 = ax.plot(time_dat[i_min:i_max],e_dat[i_min:i_max], 'k-')
    l2 = ax.plot(gap_time, inner_gap_e, 'b:')
    l3 = ax.plot(gap_time, peak_gap_e, 'b--')
    ax.set_xlabel(r'Time (years)', size =14)
    ax.set_ylabel(r'Eccentricity ($a_{KEP47}$)', size = 14)
    ax.legend([r'Small Planet', r'Inner Cavity', r'Density Peak'])
    ax.grid()
    save_path = '../Report_Plot/SmallPlanetFailur_e.png'
    repeated_plots = 0
    while (os.path.isfile(save_path)):
        save_path = '../Report_Plot/SmallPlanetFailure_e({}).png'.format(repeated_plots)
        repeated_plots += 1
    print('saving plot in {}'.format(save_path))
    fig.savefig(save_path, bbox_inches = 'tight', dpi = 320)
    plt.close(fig)

def comparefour():
    fig, ax = plt.subplots(1,1)
    lines = []
    for i in range(0,4):
        name_data = UI_helper.selectDataToPlot()
        dir = all_data_dir + name_data + '/out'
        data = pluto.Pluto(dir)
        planet = UI_helper.selectObjectToPlot(dir)
        (
            time_dat,
            a_dat,
            e_dat,
            period_dat,
            mass_dat
        ) = Data_parser_helper.getNbodyInformation_dat(name_data, planet.id)
        _ = UI_helper.selectPlottingRange(dir)
        n_min, n_max = _[0],_[1]
        t_min, t_max = nts * n_min, nts * n_max
        i_min, i_max = tools.time_split(time_dat, t_min, t_max)
        time_dat *= yrs
        t_min *= yrs
        t_max *= yrs
        lines.append(ax.plot(time_dat[i_min:i_max],a_dat[i_min:i_max]))

    lines.append(ax.plot([t_min, t_max], [kepb, kepb], 'g--'))
    lines.append(ax.plot([t_min, t_max], [a_crit, a_crit], 'r--'))
    
    ax.set_xlabel(r'Time (years)', size =14)
    ax.set_ylabel(r'Semi-Major Axis ($a_{KEP47}$)', size = 14)
    ax.legend([r'Setup 1', r'Setup 2', r'Setup 3', r'Setup 4', r'Kepler 47 b', r'Critical Radius'])
    ax.grid()
    save_path = '../Report_Plot/fourplanets.png'
    repeated_plots = 0
    while (os.path.isfile(save_path)):
        save_path = '../Report_Plot/fourplanets({}).png'.format(repeated_plots)
        repeated_plots += 1
    print('saving plot in {}'.format(save_path))
    fig.savefig(save_path, bbox_inches = 'tight', dpi = 320)
    plt.close(fig)

def disc_relaxation():
    print('Initial Conditions: \n')
    init_name = UI_helper.selectDataToPlot()
    init_dir = all_data_dir + init_name + '/out'
    init_data = pluto.Pluto(init_dir)
    N = UI_helper.selectDataFileToPlot(init_dir)
    init_var_data = init_data.primitive_variable(var, N)[0,:,:]
    init_R, init_Phi = np.meshgrid(init_data.grid['faces'][0], init_data.grid['faces'][1])
    (
        init_a,
        init_e,
        _, _, _, _, _,
    ) = gap.ellipse_find(init_R, init_Phi, init_var_data)
    print('Disc: \n')
    disc_name = UI_helper.selectDataToPlot()
    disc_dir = all_data_dir + disc_name + '/out'
    disc_data = pluto.Pluto(disc_dir)
    _ = UI_helper.selectPlottingRange(disc_dir)
    n_min, n_max = _[0],_[1]
    t_min, t_max = nts * n_min, nts * n_max

    (
        disc_gap_time,
        disc_inner_gap_e,
        disc_inner_gap_a,
        disc_peak_gap_e,
        disc_peak_gap_a,
    ) = (
        [],
        [],
        [],
        [],
        [],
    )
    for i in range(n_min, n_max+1):
        disc_gap_time.append(i*nts*yrs)
        var_data = disc_data.primitive_variable(var, i)[0,:,:]
        R, Phi = np.meshgrid(disc_data.grid['faces'][0], disc_data.grid['faces'][1])
        (
            inner_a_j,
            inner_e_j,
            _, _, _, _, _,
        ) = gap.ellipse_find(R, Phi, var_data)
        disc_inner_gap_e.append(abs(inner_e_j))
        disc_inner_gap_a.append(abs(inner_a_j))
        (
            peak_a_j,
            peak_e_j,
            _, _, _, _, _,
        ) = gap.ellipse_find_alt(R, Phi, var_data)
        disc_peak_gap_e.append(abs(peak_e_j))
        disc_peak_gap_a.append(abs(peak_a_j))

    t_min *= yrs
    t_max *= yrs

    fig, axs = plt.subplots(2,1, sharex = 'all')
    l1 = axs[0].plot(disc_gap_time, disc_inner_gap_a, 'k-')
    l2 = axs[0].plot([t_min,t_max], [kepb, kepb], 'g--')
    l3 = axs[0].plot([t_min,t_max], [init_a, init_a], 'r:')
    axs[0].set_title(r'Semi-Major Axis ($a_{KEP47}$)')

    l4 = axs[1].plot(disc_gap_time, disc_inner_gap_e, 'k-')
    l5 = axs[1].plot([t_min,t_max], [abs(init_e), abs(init_e)], 'r:')
    axs[1].set_title(r'Eccentricity')

    axs[1].set(xlabel = 'Time (years)')
    axs[0].grid()
    axs[1].grid()
    axs[0].legend([r'Cavity a', r'Kepler 47 b', r'Initial a'])
    axs[1].legend([r'Cavity e', r'Initial e'])
    save_path = '../Report_Plot/disc_relax.png'
    repeated_plots = 0
    while (os.path.isfile(save_path)):
        save_path = '../Report_Plot/disc_relax({}).png'.format(repeated_plots)
        repeated_plots += 1
    print('saving plot in {}'.format(save_path))
    fig.savefig(save_path, bbox_inches = 'tight', dpi = 320)
    plt.close(fig)
    
def oned_comparison():
    keep_plotting = True
    fig, ax = plt.subplots(1,1)
    leg = []
    norm_check = input('Normalized (y/n)? \n')
    if norm_check =='y': norm = True
    if norm_check =='n': norm = False
    while keep_plotting:
        data_name = UI_helper.selectDataToPlot()
        data_dir = all_data_dir + data_name + '/out'
        data = pluto.Pluto(data_dir)
        n = UI_helper.selectDataFileToPlot(data_dir)
        legend_name = input('Name this line: \n')
        fmt = input('fmt: \n')

        var_data = data.primitive_variable(var, n)[0,:,:]
        R, Phi = np.meshgrid(data.grid['faces'][0], data.grid['faces'][1])
        Rs = R[0,:-1]

        averaged_var = np.mean(var_data, axis=0)

        sigma_ref = data.user_def_parameters['SIGMA_REF']
        alpha_sigma = data.user_def_parameters['ALPHA_SIGMA']

        if norm: norm_var = averaged_var /(sigma_ref * np.power(Rs, -alpha_sigma))
        else: norm_var = averaged_var#/np.max(averaged_var)

        ax.plot(Rs, norm_var, fmt)
        leg.append(legend_name)

        plot_check = input('Keep Plotting (y/n)? \n')
        if plot_check == 'y': keep_plotting = True
        elif plot_check == 'n': keep_plotting = False
    
    ax.set_xlabel(r'Radius ($a_{KEPB}$)')
    ax.set_ylabel(r'Normalized surface density')
    #ax.set_xlim(0,20)
    ax.legend(leg)
    ax.grid()
    save_path = '../Report_Plot/oneddensity.png'
    repeated_plots = 0
    while (os.path.isfile(save_path)):
        save_path = '../Report_Plot/oneddensity({}).png'.format(repeated_plots)
        repeated_plots += 1
    print('saving plot in {}'.format(save_path))
    fig.savefig(save_path, bbox_inches = 'tight', dpi = 320)
    plt.close(fig)

def disc_total_mass():
    data_name = UI_helper.selectDataToPlot()
    disc_dir = all_data_dir + data_name + '/out'
    (time, m, _, _, _, _, _, _, _, _) = Data_parser_helper.get_averages_data(data_name)
    _ = UI_helper.selectPlottingRange(disc_dir)
    n_min, n_max = _[0],_[1]
    t_min, t_max = nts * n_min, nts * n_max
    i_min, i_max = tools.time_split(time, t_min, t_max)
    time *= yrs
    t_min *= yrs
    t_max *= yrs
    plt.plot(time[i_min:i_max], m[i_min:i_max])
    plt.grid()
    plt.ylabel(r'Total Disc mass ($M_{KEP47}$)')
    plt.xlabel(r'Time (years)')
    save_path = '../Report_Plot/disc_mass.png'
    repeated_plots = 0
    while (os.path.isfile(save_path)):
        save_path = '../Report_Plot/disc_mass({}).png'.format(repeated_plots)
        repeated_plots += 1
    print('saving plot in {}'.format(save_path))
    plt.savefig(save_path, bbox_inches = 'tight', dpi = 320)
    plt.close()

def two_planets():
    data_name = UI_helper.selectDataToPlot()
    data_dir = all_data_dir + data_name + '/out'
    data = pluto.Pluto(data_dir)
    p1 = UI_helper.selectObjectToPlot(data_dir)
    (
        time_dat_1,
        a_dat_1,
        e_dat_1,
        period_dat_1,
        mass_dat_1,
    ) = Data_parser_helper.getNbodyInformation_dat(data_name, p1.id)
    p2 = UI_helper.selectObjectToPlot(data_dir)
    (
        time_dat_2,
        a_dat_2,
        e_dat_2,
        period_dat_2,
        mass_dat_2,
    ) = Data_parser_helper.getNbodyInformation_dat(data_name, p2.id)
    _ = UI_helper.selectPlottingRange(data_dir)
    n_min, n_max = _[0],_[1]
    t_min, t_max = nts * n_min, nts * n_max
    i_min_1, i_max_1 = tools.time_split(time_dat_1, t_min, t_max)
    i_min_2, i_max_2 = tools.time_split(time_dat_2, t_min, t_max)
    time_dat_1 *= yrs
    time_dat_2 *= yrs
    t_min *= yrs
    t_max *= yrs
    (
        gap_time,
        inner_gap_e,
        inner_gap_a,
        peak_gap_e,
        peak_gap_a,
    ) = (
        [],
        [],
        [],
        [],
        [],
    )
    for i in range(n_min, n_max+1):
        gap_time.append(i*nts*yrs)
        var_data = data.primitive_variable(var, i)[0,:,:]
        R, Phi = np.meshgrid(data.grid['faces'][0], data.grid['faces'][1])
        (
            inner_a_j,
            inner_e_j,
            _, _, _, _, _,
        ) = gap.ellipse_find(R, Phi, var_data)
        inner_gap_e.append(abs(inner_e_j))
        inner_gap_a.append(abs(inner_a_j))
        (
            peak_a_j,
            peak_e_j,
            _, _, _, _, _,
        ) = gap.ellipse_find_alt(R, Phi, var_data)
        peak_gap_e.append(abs(peak_e_j))
        peak_gap_a.append(abs(peak_a_j))
    fig, ax = plt.subplots(1,1)
    l1 = ax.plot(time_dat_1[i_min_1:i_max_1], a_dat_1[i_min_1:i_max_1], 'k-')
    l2 = ax.plot(time_dat_2[i_min_2:i_max_2], a_dat_2[i_min_2:i_max_2], 'm-')
    l3 = ax.plot(gap_time, inner_gap_a, 'b:')
    l4 = ax.plot(gap_time, peak_gap_a, 'b--')
    l5 = ax.plot([t_min, t_max], [kepb, kepb], 'g--')
    l6 = ax.plot([t_min, t_max], [a_crit, a_crit], 'r--')
    ax.grid()
    ax.set_ylabel(r'Semi-Major Axis ($a_{KEP47}$)')
    ax.set_xlabel(r'Time (years)')
    ax.legend(['Large Planet', 'Small Planet', 'Inner Cavity', 'Disc Peak', 'Kepler 47 b', 'Critical Radius'])
    save_path = '../Report_Plot/two_planets.png'
    repeated_plots = 0
    while (os.path.isfile(save_path)):
        save_path = '../Report_Plot/two_planets({}).png'.format(repeated_plots)
        repeated_plots += 1
    print('saving plot in {}'.format(save_path))
    plt.savefig(save_path, bbox_inches = 'tight', dpi = 320)
    plt.close()

def plot_resonance_1():
    data_name = UI_helper.selectDataToPlot()
    data_dir = all_data_dir + data_name + '/out'
    data = pluto.Pluto(data_dir)
    p1 = UI_helper.selectObjectToPlot(data_dir)
    (
        time_dat_1,
        a_dat_1,
        e_dat_1,
        omega_dat_1,
        anomoly_dat_1,
    ) = Data_parser_helper.getNbodyInformation_out(data_name, p1)
    p2 = UI_helper.selectObjectToPlot(data_dir)
    (
        time_dat_2,
        a_dat_2,
        e_dat_2,
        omega_dat_2,
        anomoly_dat_2,
    ) = Data_parser_helper.getNbodyInformation_out(data_name, p2)
    _ = UI_helper.selectPlottingRange(data_dir)
    n_min, n_max = _[0],_[1]
    t_min, t_max = nts * n_min, nts * n_max
    # i_min_1, i_max_1 = tools.time_split(time_dat_1, t_min, t_max)
    # i_min_2, i_max_2 = tools.time_split(time_dat_2, t_min, t_max)
    time_dat_1 *= yrs
    time_dat_2 *= yrs
    t_min *= yrs
    t_max *= yrs
    (
        gap_time,
        inner_gap_e,
        inner_gap_a,
        peak_gap_e,
        peak_gap_a,
    ) = (
        [],
        [],
        [],
        [],
        [],
    )
    for i in range(n_min, n_max+1):
        gap_time.append(i*nts*yrs)
        var_data = data.primitive_variable(var, i)[0,:,:]
        R, Phi = np.meshgrid(data.grid['faces'][0], data.grid['faces'][1])
        (
            inner_a_j,
            inner_e_j,
            _, _, _, _, _,
        ) = gap.ellipse_find(R, Phi, var_data)
        inner_gap_e.append(abs(inner_e_j))
        inner_gap_a.append(abs(inner_a_j))
        (
            peak_a_j,
            peak_e_j,
            _, _, _, _, _,
        ) = gap.ellipse_find_alt(R, Phi, var_data)
        peak_gap_e.append(abs(peak_e_j))
        peak_gap_a.append(abs(peak_a_j))
    
    a_1 = a_dat_1[n_min:n_max+1]
    a_2 = a_dat_2[n_min:n_max+1]

    r11_a = np.power(a_1/a_1, 3/2)
    r12_a = np.power(a_2/a_1, 3/2)
    r1i_a = np.power(inner_gap_a/a_1, 3/2)
    r1p_a = np.power(peak_gap_a/a_1, 3/2)

    fig, ax = plt.subplots(1,1)
    l1 = ax.plot(gap_time, r11_a, 'k-')
    l2 = ax.plot(gap_time, r12_a, 'm-')
    l3 = ax.plot(gap_time, r1i_a, 'b:')
    l4 = ax.plot(gap_time, r1p_a, 'b--')
    l5 = ax.plot([gap_time[n_min], gap_time[n_max]], [4/3, 4/3], 'g--')
    l6 = ax.plot([gap_time[n_min], gap_time[n_max]], [8/3, 8/3], 'g--')
    ax.grid()
    ax.set_ylabel(r'Ratio of orbital period')
    ax.set_xlabel(r'Time (years)')
    ax.legend(['Large Planet', 'Small Planet', 'Inner Cavity', 'Disc Peak'])
    save_path = '../Report_Plot/resonance.png'
    repeated_plots = 0
    while (os.path.isfile(save_path)):
        save_path = '../Report_Plot/resonance({}).png'.format(repeated_plots)
        repeated_plots += 1
    print('saving plot in {}'.format(save_path))
    plt.savefig(save_path, bbox_inches = 'tight', dpi = 320)
    plt.close()

def plot_resonance_2():
    data_name = UI_helper.selectDataToPlot()
    data_dir = all_data_dir + data_name + '/out'
    data = pluto.Pluto(data_dir)
    p1 = UI_helper.selectObjectToPlot(data_dir)
    (
        time_dat_1,
        a_dat_1,
        e_dat_1,
        omega_dat_1,
        anomoly_dat_1,
    ) = Data_parser_helper.getNbodyInformation_out(data_name, p1)
    p2 = UI_helper.selectObjectToPlot(data_dir)
    (
        time_dat_2,
        a_dat_2,
        e_dat_2,
        omega_dat_2,
        anomoly_dat_2,
    ) = Data_parser_helper.getNbodyInformation_out(data_name, p2)
    _ = UI_helper.selectPlottingRange(data_dir)
    n_min, n_max = _[0],_[1]
    t_min, t_max = nts * n_min, nts * n_max
    # i_min_1, i_max_1 = tools.time_split(time_dat_1, t_min, t_max)
    # i_min_2, i_max_2 = tools.time_split(time_dat_2, t_min, t_max)
    time_dat_1 *= yrs
    time_dat_2 *= yrs
    t_min *= yrs
    t_max *= yrs
    (
        gap_time,
        inner_gap_e,
        inner_gap_a,
        peak_gap_e,
        peak_gap_a,
    ) = (
        [],
        [],
        [],
        [],
        [],
    )
    for i in range(n_min, n_max+1):
        gap_time.append(i*nts*yrs)
        var_data = data.primitive_variable(var, i)[0,:,:]
        R, Phi = np.meshgrid(data.grid['faces'][0], data.grid['faces'][1])
        (
            inner_a_j,
            inner_e_j,
            _, _, _, _, _,
        ) = gap.ellipse_find(R, Phi, var_data)
        inner_gap_e.append(abs(inner_e_j))
        inner_gap_a.append(abs(inner_a_j))
        (
            peak_a_j,
            peak_e_j,
            _, _, _, _, _,
        ) = gap.ellipse_find_alt(R, Phi, var_data)
        peak_gap_e.append(abs(peak_e_j))
        peak_gap_a.append(abs(peak_a_j))
    
    a_1 = a_dat_1[n_min:n_max+1]
    a_2 = a_dat_2[n_min:n_max+1]

    r11_a = np.power(a_1/a_1, 3/2)
    r22_a = np.power(a_2/a_2, 3/2)
    r12_a = np.power(a_2/a_1, 3/2)
    r1i_a = np.power(inner_gap_a/a_1, 3/2)
    r1p_a = np.power(peak_gap_a/a_1, 3/2)
    ri2_a = np.power(inner_gap_a/a_2, -3/2)

    fig, ax = plt.subplots(1,1)
    l1 = ax.plot(gap_time, r22_a, 'k-')
    #l2 = ax.plot(gap_time, r12_a, 'm-')
    l2 = ax.plot(gap_time, ri2_a, 'b:')
    #l4 = ax.plot(gap_time, r1p_a, 'b--')
    #l5 = ax.plot([gap_time[n_min], gap_time[n_max]], [4/3, 4/3], 'g--')
    l6 = ax.plot([gap_time[n_min], gap_time[n_max]], [2, 2], 'g--')
    ax.grid()
    ax.set_ylabel(r'Ratio of orbital period compared to inner cavity')
    ax.set_xlabel(r'Time (years)')
    ax.legend(['Small Planet', 'Inner Cavity'])
    save_path = '../Report_Plot/resonance.png'
    repeated_plots = 0
    while (os.path.isfile(save_path)):
        save_path = '../Report_Plot/resonance({}).png'.format(repeated_plots)
        repeated_plots += 1
    print('saving plot in {}'.format(save_path))
    plt.savefig(save_path, bbox_inches = 'tight', dpi = 320)
    plt.close()

def plot_2_planets():
    name1 = UI_helper.selectDataToPlot()
    dir1 = all_data_dir + name1 + '/out'
    data1 = pluto.Pluto(dir1)
    p1 = UI_helper.selectObjectToPlot(dir1)
    n1_min, n1_max = UI_helper.selectPlottingRange(dir1)
    name2 = UI_helper.selectDataToPlot()
    dir2 = all_data_dir + name2 + '/out'
    data2 = pluto.Pluto(dir2)
    p2 = UI_helper.selectObjectToPlot(dir2)
    p3 = UI_helper.selectObjectToPlot(dir2)
    n2_min, n2_max = UI_helper.selectPlottingRange(dir2)
    (
        time_dat_1,
        a_dat_1,
        e_dat_1,
        omega_dat_1,
        anomoly_dat_1,
    ) = Data_parser_helper.getNbodyInformation_dat(name1, p1.id)
    (
        time_dat_2,
        a_dat_2,
        e_dat_2,
        omega_dat_2,
        anomoly_dat_2,
    ) = Data_parser_helper.getNbodyInformation_dat(name2, p2.id)
    (
        time_dat_3,
        a_dat_3,
        e_dat_3,
        omega_dat_3,
        anomoly_dat_3,
    ) = Data_parser_helper.getNbodyInformation_dat(name2, p3.id)
    t1_min, t1_max = nts * n1_min, nts * n1_max
    t2_min, t2_max = nts * n2_min, nts * n2_max
    i1_min, i1_max = tools.time_split(time_dat_1, t1_min, t1_max)
    i2_min, i2_max = tools.time_split(time_dat_2, t2_min, t2_max)
    i3_min, i3_max = tools.time_split(time_dat_3, t2_min, t2_max)
    time_dat_1 *= yrs
    time_dat_2 *= yrs
    time_dat_3 *= yrs
    t1_min *= yrs
    t1_max *= yrs
    t2_min *= yrs
    t2_max *= yrs
    (
        gap_time1,
        gap_time2,
        inner1_gap_e,
        inner1_gap_a,
        inner2_gap_e,
        inner2_gap_a,
    ) = (
        [],
        [],
        [],
        [],
        [],
        [],
    )
    for i in range(n1_min, n1_max+1):
        gap_time1.append(i*nts*yrs)
        var_data1 = data1.primitive_variable(var, i)[0,:,:]
        R, Phi = np.meshgrid(data1.grid['faces'][0], data1.grid['faces'][1])
        (
            inner1_a_j,
            inner1_e_j,
            _, _, _, _, _,
        ) = gap.ellipse_find(R, Phi, var_data1)
        inner1_gap_e.append(abs(inner1_e_j))
        inner1_gap_a.append(abs(inner1_a_j))
    
    for i in range(n2_min, n2_max+1):
        gap_time2.append(i*nts*yrs)
        var_data2 = data2.primitive_variable(var, i)[0,:,:]
        R, Phi = np.meshgrid(data2.grid['faces'][0], data2.grid['faces'][1])
        (
            inner2_a_j,
            inner2_e_j,
            _, _, _, _, _,
        ) = gap.ellipse_find(R, Phi, var_data2)
        inner2_gap_e.append(abs(inner2_e_j))
        inner2_gap_a.append(abs(inner2_a_j))

    offset = int(input('Offset? \n'))
    offset_t = offset * nts * yrs
    
    fig, ax = plt.subplots(1,1)
    l1 = ax.plot(time_dat_1[i1_min:i1_max], a_dat_1[i1_min:i1_max], 'k-')
    l2 = ax.plot(time_dat_2[i2_min:i2_max] + offset_t, a_dat_2[i2_min:i2_max], 'g-')
    l3 = ax.plot(time_dat_3[i2_min:i2_max] + offset_t, a_dat_3[i2_min:i2_max], 'm-')
    l4 = ax.plot(gap_time1, inner1_gap_a, 'k:')
    l5 = ax.plot(np.array(gap_time2) + offset_t, inner2_gap_a, 'g:')
    l6 = ax.plot([t1_min,t1_max], [a_crit, a_crit], 'r--')
    ax.grid()
    ax.set_ylabel(r'Semi-Major Axis ($a_{KEP47}$)')
    ax.set_xlabel(r'Time (years)')
    ax.legend(['Single planet', 'Planet with smaller companion', 'Smaller Companion', 'Single planet cavity', 'Accompanied planet cavity', 'Critical Radius'])
    save_path = '../Report_Plot/withwithout.png'
    repeated_plots = 0
    while (os.path.isfile(save_path)):
        save_path = '../Report_Plot/withwithout({}).png'.format(repeated_plots)
        repeated_plots += 1
    print('saving plot in {}'.format(save_path))
    plt.savefig(save_path, bbox_inches = 'tight', dpi = 320)
    plt.close()

def plot_planet_mass():
    

    name_data = UI_helper.selectDataToPlot()
    dir = all_data_dir + name_data + '/out'
    data = pluto.Pluto(dir)
    planet = UI_helper.selectObjectToPlot(dir)
    (
        time_dat,
        a_dat,
        e_dat,
        period_dat,
        mass_dat
    ) = Data_parser_helper.getNbodyInformation_dat(name_data, planet.id)
    _ = UI_helper.selectPlottingRange(dir)
    n_min, n_max = _[0],_[1]
    nts = int(input('nts \n'))
    (time, m, _, _, _, _, _, _, _, _) = Data_parser_helper.get_averages_data(name_data)
    t_min, t_max = nts * n_min, nts * n_max
    i_min, i_max = tools.time_split(time_dat, t_min, t_max)
    id_min, id_max = tools.time_split(time, t_min, t_max)

    # (
    #     gap_time,
    #     inner_gap_e,
    #     inner_gap_a,
    #     peak_gap_e,
    #     peak_gap_a,
    # ) = (
    #     [],
    #     [],
    #     [],
    #     [],
    #     [],
    # )
    # for i in range(n_min, n_max+1):
    #     gap_time.append(i*nts*yrs)
    #     var_data = data.primitive_variable(var, i)[0,:,:]
    #     R, Phi = np.meshgrid(data.grid['faces'][0], data.grid['faces'][1])
    #     (
    #         inner_a_j,
    #         inner_e_j,
    #         _, _, _, _, _,
    #     ) = gap.ellipse_find(R, Phi, var_data)
    #     inner_gap_e.append(abs(inner_e_j))
    #     inner_gap_a.append(abs(inner_a_j))
    #     (
    #         peak_a_j,
    #         peak_e_j,
    #         _, _, _, _, _,
    #     ) = gap.ellipse_find_alt(R, Phi, var_data)
    #     peak_gap_e.append(abs(peak_e_j))
    #     peak_gap_a.append(abs(peak_a_j))

    time_dat *= yrs
    t_min *= yrs
    t_max *= yrs

    fig, ax = plt.subplots(1,1)
    plt.plot(time_dat[i_min:i_max], mass_dat[i_min:i_max]- 0.001, 'r-')
    plt.plot(time[id_min:id_max]*yrs, (m[0] - m[id_min:id_max]), 'r--')
    #plt.plot(time[id_min:id_max]*yrs, ((mass_dat-0.001)/(m[0]- m))[id_min:id_max])
    plt.ylabel(r'Accreted mass ($M_{KEP47}$)')
    plt.xlabel(r'Time (years)')
    plt.legend([r'Mass gained by planet', r'Mass lost by Disc'])
    plt.text(600,
             0.0,
             r'$\tau_{ACC} = 30$ years',
             color = 'k',
             size =20)
    plt.grid()

    save_path = '../Report_Plot/accretion.png'
    repeated_plots = 0
    while (os.path.isfile(save_path)):
        save_path = '../Report_Plot/accretion({}).png'.format(repeated_plots)
        repeated_plots += 1
    print('saving plot in {}'.format(save_path))
    plt.savefig(save_path, bbox_inches = 'tight', dpi = 320)
    plt.close()
    print(len(time))
    print(len(time_dat))

def compare_accretions():
    keep_plotting = True
    fig, ax = plt.subplots(1,1)
    leg = []
    while keep_plotting:
        name_data = UI_helper.selectDataToPlot()
        dir = all_data_dir + name_data + '/out'
        data = pluto.Pluto(dir)
        planet = UI_helper.selectObjectToPlot(dir)
        (
            time_dat,
            a_dat,
            e_dat,
            period_dat,
            mass_dat
        ) = Data_parser_helper.getNbodyInformation_dat(name_data, planet.id)
        _ = UI_helper.selectPlottingRange(dir)
        n_min, n_max = _[0],_[1]
        legend_name = str(input('Name this line: \n'))
        leg.append(legend_name)
        fmt = str(input('fmt: \n'))
        nts = int(input('nts \n'))
        (time, m, _, _, _, _, _, _, _, _) = Data_parser_helper.get_averages_data(name_data)
        t_min, t_max = nts * n_min, nts * n_max
        i_min, i_max = tools.time_split(time_dat, t_min, t_max)
        ax.plot(time[i_min:i_max]*yrs, ((mass_dat-0.001)/(m[0]- m))[i_min:i_max], fmt)
        plot_check = input('Keep Plotting (y/n)? \n')
        if plot_check == 'y': keep_plotting = True
        elif plot_check == 'n': keep_plotting = False
    
    ax.set_ylabel(r'Accretion$_{planet}$/accretion$_{disc}$')
    ax.set_xlabel(r'Time (years)')
    ax.legend(leg)
    ax.grid()

    save_path = '../Report_Plot/accretion_comparison.png'
    repeated_plots = 0
    while (os.path.isfile(save_path)):
        save_path = '../Report_Plot/accretion_comparison({}).png'.format(repeated_plots)
        repeated_plots += 1
    print('saving plot in {}'.format(save_path))
    plt.savefig(save_path, bbox_inches = 'tight', dpi = 320)
    plt.close()
    print(len(time))
    print(len(time_dat))

def several_planets():
    keep_plotting = True
    fig, ax = plt.subplots(1,1)
    leg = []
    n_maxes = []
    while keep_plotting:
        name_data = UI_helper.selectDataToPlot()
        dir = all_data_dir + name_data + '/out'
        planet = UI_helper.selectObjectToPlot(dir)
        (
            time_dat,
            a_dat,
            e_dat,
            period_dat,
            mass_dat
        ) = Data_parser_helper.getNbodyInformation_dat(name_data, planet.id)
        _ = UI_helper.selectPlottingRange(dir)
        n_min, n_max = _[0],_[1]
        legend_name = str(input('Name this line: \n'))
        leg.append(legend_name)
        fmt = str(input('fmt: \n'))
        nts = int(input('nts \n'))
        t_min, t_max = nts * n_min, nts * n_max
        i_min, i_max = tools.time_split(time_dat, t_min, t_max)
        ax.plot(time_dat[i_min:i_max]*yrs, a_dat[i_min:i_max], fmt)
        plot_check = input('Keep Plotting (y/n)? \n')
        n_maxes.append(n_max)
        if plot_check == 'y': keep_plotting = True
        elif plot_check == 'n': keep_plotting = False
    ax.grid()
    n_max = np.max(n_maxes)
    n_min = 0
    n_max *= nts * yrs
    ax.plot([n_min,n_max], [kepb, kepb])
    leg.append('Kepler 47 b')
    ax.legend(leg)
    ax.set_ylabel(r'Semi-major axis ($a_{KEP47}$)')
    ax.set_xlabel(r'Time (years)')

if __name__ == '__main__':
    plotters = [
        small_planet_typeI_blocked,
        comparefour,
        disc_relaxation,
        oned_comparison,
        disc_total_mass,
        two_planets,
        plot_resonance_1,
        plot_resonance_2,
        plot_2_planets,
        plot_planet_mass,
        compare_accretions,
        several_planets,
        ]
    func_index = UI_helper.selectFunctionsToRun(plotters)
    eval('{}()'.format(plotters[func_index].__name__))
