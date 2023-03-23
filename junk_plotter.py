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
from celluloid import Camera

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

plots_dir = '../Poster_plots/'
n = 200
n_min = 0
n_max = n
t_min = nts * n_min
t_max = nts * n_max

a_crit = 2.31
kepb = 3.53

a3h3_name = 'acc_alpha-3h.03_MJ1.0_aB13.0_2/'
a4h4_name = 'acc_47a-4_MJ1.0_aB13.0_2/'
a3h4_name = 'acc_MJ1.0_aB13.0_Ve-3_2/'
a4h3_name = 'acc_MJ1.0_aB13.0_Ve-4_2/'
a3h4_na_name = 'MJ1.0_aB13.0_Ve-3_restart/'

a3h3_dir = all_data_dir + a3h3_name + 'out/'
a3h4_dir = all_data_dir + a3h4_name + 'out/'
a4h3_dir = all_data_dir + a4h3_name + 'out/'
a4h4_dir = all_data_dir + a4h4_name + 'out/'
a3h4_na_dir = all_data_dir + a3h4_na_name + 'out/'

a3h3_data = pluto.Pluto(a3h3_dir)
a3h4_data = pluto.Pluto(a3h4_dir)
a4h3_data = pluto.Pluto(a4h3_dir)
a4h4_data = pluto.Pluto(a4h4_dir)
a3h4_na_data = pluto.Pluto(a3h4_na_dir)

(
    a3h3_time_dat,
    a3h3_a_dat,
    a3h3_e_dat,
    a3h3_period_dat,
    a3h3_mass_dat
) = Data_parser_helper.getNbodyInformation_dat(a3h3_name, 2)
(
    a3h4_time_dat,
    a3h4_a_dat,
    a3h4_e_dat,
    a3h4_period_dat,
    a3h4_mass_dat
) = Data_parser_helper.getNbodyInformation_dat(a3h4_name, 2)
(
    a4h3_time_dat,
    a4h3_a_dat,
    a4h3_e_dat,
    a4h3_period_dat,
    a4h3_mass_dat
) = Data_parser_helper.getNbodyInformation_dat(a4h3_name, 2)
(
    a4h4_time_dat,
    a4h4_a_dat,
    a4h4_e_dat,
    a4h4_period_dat,
    a4h4_mass_dat
) = Data_parser_helper.getNbodyInformation_dat(a4h4_name, 2)
(
    a3h4_na_time_dat,
    a3h4_na_a_dat,
    a3h4_na_e_dat,
    a3h4_na_period_dat,
    a3h4_na_mass_dat
) = Data_parser_helper.getNbodyInformation_dat(a3h4_na_name, 2)

#a3h3_Nchar_dat = Data_parser_helper.getNbodyInformation_dat(a3h3_name, 2)
#a3h4_Nchar_dat = Data_parser_helper.getNbodyInformation_dat(a3h4_name, 2)
#a4h3_Nchar_dat = Data_parser_helper.getNbodyInformation_dat(a4h3_name, 2)
#a4h4_Nchar_dat = Data_parser_helper.getNbodyInformation_dat(a4h4_name, 2)

def plot_141122_vs():
    name = '141122_1_restart/'
    dir = all_data_dir + name + 'out/'
    data = pluto.Pluto(dir)
    (
        time_dat,
        a_dat,
        e_dat,
        period_dat,
        mass_dat
    ) = Data_parser_helper.getNbodyInformation_dat(name, 2)
    i_min, i_max = tools.time_split(time_dat, t_min, t_max)
    
    gap_time = []
    gap_e = []
    gap_a = []
    for j in range(n_min, n_max+1):
        gap_time.append(j*nts)
        var_data = data.primitive_variable(var, j)[0,:,:]
        R, Phi = np.meshgrid(data.grid['faces'][0], data.grid['faces'][1])
        (
            a_j,
            e_j,
            _, _, _, _, _
        ) = gap.ellipse_find(R, Phi, var_data)
        gap_e.append(abs(e_j))
        gap_a.append(abs(a_j))

    fig, ax1 = plt.subplots(1,1)
#    ax1.grid(color = '#AAAAFF')
    l1, = ax1.plot(time_dat[i_min:i_max],a_dat[i_min:i_max], 'b')
    ax2 = ax1.twinx()
    plotp.poster()
#    ax2.grid(color = '#FFAAAA')
    ax1.set_ylim(3.0,13.0)
    ax1.set_ylabel(r'Semi-Major Axis ($a_{bin}$)', size = 14, color = 'b')
    ax1.set_xlabel(r'Time ($T_{bin}$)', size = 14)
    ax2.set_ylabel(r'Eccentricity ($e$)', size = 14, color = 'r')
    l2, = ax2.plot(time_dat[i_min:i_max],e_dat[i_min:i_max], 'r')
    l3, = ax1.plot(gap_time, gap_a, 'b:')
    l4, = ax2.plot(gap_time, gap_e, 'r:')
    l5, = ax1.plot([t_min,t_max], [kepb,kepb], 'g--', linewidth = 2)
    ax2.legend([l1,l2,l3,l4,l5], ['planet a', 'planet e','gap a', 'gap e', 'Kep47b a'])
    save_path = '../Poster_plots/Finished_plot_1.png'
    repeated_plots = 0
    while (os.path.isfile(save_path)):
        save_path = '../Poster_plots/Finished_plot_1({}).png'.format(repeated_plots)
        repeated_plots += 1
    
    print('saving plot in {}'.format(save_path))
    fig.savefig(save_path, bbox_inches = 'tight', dpi = 320)
    plt.close(fig)

def plot_4_N_dat():
    a3h3_i_min, a3h3_i_max = tools.time_split(a3h3_time_dat, t_min, t_max)
    a4h3_i_min, a4h3_i_max = tools.time_split(a4h3_time_dat, t_min, t_max)
    a3h4_i_min, a3h4_i_max = tools.time_split(a3h4_time_dat, t_min, t_max)
    a4h4_i_min, a4h4_i_max = tools.time_split(a4h4_time_dat, t_min, t_max)
    a3h4_na_i_min, a3h4_na_i_max = tools.time_split(a3h4_na_time_dat, t_min, t_max)
    
    fig = plt.figure()
    plotp.poster()
#   plt.title(r'Planet Migration ($M_P = M_J$)')
    plt.plot(a4h3_time_dat[a4h3_i_min:a4h3_i_max], a4h3_a_dat[a4h3_i_min:a4h3_i_max], label = r'setup 1')
    plt.plot(a3h3_time_dat[a3h3_i_min:a3h3_i_max], a3h3_a_dat[a3h3_i_min:a3h3_i_max], label = r'setup 2')
    plt.plot(a4h4_time_dat[a4h4_i_min:a4h4_i_max], a4h4_a_dat[a4h4_i_min:a4h4_i_max], label = r'setup 3')
    plt.plot(a3h4_time_dat[a3h4_i_min:a3h4_i_max], a3h4_a_dat[a3h4_i_min:a3h4_i_max], label = r'setup 4')
    plt.plot(a3h4_na_time_dat[a3h4_na_i_min:a3h4_na_i_max], a3h4_na_a_dat[a3h4_na_i_min:a3h4_na_i_max], label = r'setup 4 no accretion')
    plt.xlabel('Time ($T_{bin}$)', size = 14)
    plt.xticks(np.arange(0,10001, step = 2500))
    plt.ylabel('Semi-Major Axis ($a_{bin}$)', size = 14)
    #plt.grid()
    plt.legend()
    fig.tight_layout
    save_path = '../Poster_plots/Finished_plot_2.png'
    repeated_plots = 0
    while (os.path.isfile(save_path)):
        save_path = '../Poster_plots/Finished_plot_2({}).png'.format(repeated_plots)
        repeated_plots += 1
    
    print('Saving plot in {}'.format(save_path))
    fig.savefig(save_path,bbox_inches = 'tight', dpi = 320)
    plt.close(fig)

def plot_4_2ds_kms():
    name = 'MJ1.0_aB13.0_Ve-4_restart/'
    dir = all_data_dir + name + 'out/'
    #directories = Navigation_helper.Directories(name)
    data = pluto.Pluto(dir)

    var_data_0 = data.primitive_variable(var, 0)
    var_data_1 = data.primitive_variable(var, 50)
    var_data_2 = data.primitive_variable(var, 100)
    var_data_3 = data.primitive_variable(var, 150)
    r, phi = data.grid['centers'].X1, data.grid['centers'].X2

    fig, axs = plt.subplots(2,2)

    axs[0].set_ylabel(r'$y\;[a_\mathrm{b}]$',fontsize=fs)
    axs[0].set_xlabel(r'$x\;[a_\mathrm{b}]$',fontsize=fs)
    axs[1].set_ylabel(r'$y\;[a_\mathrm{b}]$',fontsize=fs)
    axs[1].set_xlabel(r'$x\;[a_\mathrm{b}]$',fontsize=fs)
    axs[2].set_ylabel(r'$y\;[a_\mathrm{b}]$',fontsize=fs)
    axs[2].set_xlabel(r'$x\;[a_\mathrm{b}]$',fontsize=fs)
    axs[3].set_ylabel(r'$y\;[a_\mathrm{b}]$',fontsize=fs)
    axs[3].set_xlabel(r'$x\;[a_\mathrm{b}]$',fontsize=fs)

    r_cylinder = data.grid['faces'][0][0]
    cylinder = patches.Arc((0,0), width = 2*r_cylinder,
                           height=2*r_cylinder)
    axs[0].add_patch(cylinder)
    axs[1].add_patch(cylinder)
    axs[2].add_patch(cylinder)
    axs[3].add_patch(cylinder)

    R, Phi = np.meshgrid(data.grid['faces'][0], data.grid['faces'][1])

    X = np.multiply(R,np.cos(Phi))
    Y = np.multiply(R,np.sin(Phi))

    a_0, num_e_0, phi_0_0, x_0_0, y_0_0, x_ell_0, y_ell_0 = gap.ellipse_find(R, Phi, var_data_0)
    a_1, num_e_1, phi_0_1, x_0_1, y_0_1, x_ell_1, y_ell_1 = gap.ellipse_find(R, Phi, var_data_1)
    a_2, num_e_2, phi_0_2, x_0_2, y_0_2, x_ell_2, y_ell_2 = gap.ellipse_find(R, Phi, var_data_2)
    a_3, num_e_3, phi_0_3, x_0_3, y_0_3, x_ell_3, y_ell_3 = gap.ellipse_find(R, Phi, var_data_3)

    corr = np.sqrt(1+(-1.0-0.5)*0.05*0.05)
    plot_0 = axs[0].pcolor(X*a_bin, Y*a_bin, np.log10(var_data_0/np.max(var_data_0),cmap = 'gist_heat, vmin = -2'))
    plot_1 = axs[1].pcolor(X*a_bin, Y*a_bin, np.log10(var_data_1/np.max(var_data_1),cmap = 'gist_heat, vmin = -2'))
    plot_2 = axs[2].pcolor(X*a_bin, Y*a_bin, np.log10(var_data_2/np.max(var_data_2),cmap = 'gist_heat, vmin = -2'))
    plot_3 = axs[3].pcolor(X*a_bin, Y*a_bin, np.log10(var_data_3/np.max(var_data_3),cmap = 'gist_heat, vmin = -2'))
    colourbar_label = r'$log\Sigma\;\left[\mathrm{g}\,/\mathrm{cm}^{-2}\right]$'

    cb_3 = plt.colorbar(plot_3, orientation='vertical')
    cb_3.set_label(colourbar_label, fontsize = fs)

    axs[0].plot(0,0, 'k+', ms = 3)
    axs[1].plot(0,0, 'k+', ms = 3)
    axs[2].plot(0,0, 'k+', ms = 3)
    axs[3].plot(0,0, 'k+', ms = 3)

    time_text_0 = axs[0].text(0.95, 0.95, r't={} yr'.format(0*nts),
                              color='b',
                              size='x-large',
                              ha='right',
                              va='top',
                              transform=axs[0].transAxes)
    time_text_1 = axs[1].text(0.95, 0.95, r't={} yr'.format(50*nts),
                              color='b',
                              size='x-large',
                              ha='right',
                              va='top',
                              transform=axs[1].transAxes)
    time_text_2 = axs[2].text(0.95, 0.95, r't={} yr'.format(100*nts),
                              color='b',
                              size='x-large',
                              ha='right',
                              va='top',
                              transform=axs[2].transAxes)
    time_text_3 = axs[3].text(0.95, 0.95, r't={} yr'.format(150*nts),
                              color='b',
                              size='x-large',
                              ha='right',
                              va='top',
                              transform=axs[3].transAxes)
    
    axs[0].plot(x_ell_0*a_bin, y_ell_0*a_bin, 'w--', lw=1.5)
    axs[1].plot(x_ell_1*a_bin, y_ell_1*a_bin, 'w--', lw=1.5)
    axs[2].plot(x_ell_2*a_bin, y_ell_2*a_bin, 'w--', lw=1.5)
    axs[3].plot(x_ell_3*a_bin, y_ell_3*a_bin, 'w--', lw=1.5)

    binary_text = axs[0].text(0.05, 0.15,
            r'$e_\mathrm{{gap}} = {0:.2f}$'.format(num_e_0)+'\n'
            +r'$a_\mathrm{{gap}} = {0:.2f}\,a_\mathrm{{b}}$'.format(a_0),
            color='w',
            size='x-large',
            ha='left',
            va='top',
            transform=axs[0].transAxes)
    binary_text = axs[1].text(0.05, 0.15,
            r'$e_\mathrm{{gap}} = {0:.2f}$'.format(num_e_1)+'\n'
            +r'$a_\mathrm{{gap}} = {0:.2f}\,a_\mathrm{{b}}$'.format(a_1),
            color='w',
            size='x-large',
            ha='left',
            va='top',
            transform=axs[1].transAxes)
    binary_text = axs[2].text(0.05, 0.15,
            r'$e_\mathrm{{gap}} = {0:.2f}$'.format(num_e_2)+'\n'
            +r'$a_\mathrm{{gap}} = {0:.2f}\,a_\mathrm{{b}}$'.format(a_2),
            color='w',
            size='x-large',
            ha='left',
            va='top',
            transform=axs[2].transAxes)
    binary_text = axs[3].text(0.05, 0.15,
            r'$e_\mathrm{{gap}} = {0:.2f}$'.format(num_e_3)+'\n'
            +r'$a_\mathrm{{gap}} = {0:.2f}\,a_\mathrm{{b}}$'.format(a_3),
            color='w',
            size='x-large',
            ha='left',
            va='top',
            transform=axs[3].transAxes)
    fig.tight_layout()
    plotp.poster()
    save_path = '../Poster_plots/junk/test_plot.png'
    repeated_plots = 0
    while(os.path.isfile(save_path)):
        save_path = '../Poster_plots/junk/test_plot({})'.format(repeated_plots)
        repeated_plots+=1
    
    print('Saving plot in {}'.format(save_path))
    plt.savefig(save_path)
    fig.close()

def plot_4_2ds():
    name = 'MJ1.0_aB13.0_Ve-3_restart/'
    dir = all_data_dir + name + 'out/'
    #directories = Navigation_helper.Directories(name)
    data = pluto.Pluto(dir)
    fig, axs = plt.subplots(2,2)
    TwoD_sigma.plot_the_data(0, dir, data, axs[0,0])
    TwoD_sigma.plot_the_data(50, dir, data, axs[0,1])
    TwoD_sigma.plot_the_data(100, dir, data, axs[1,0])
    TwoD_sigma.plot_the_data(150, dir, data, axs[1,1])

    fig.tight_layout()
    plotp.poster()
    save_path = '../Poster_plots/junk/test_plot.png'
    repeated_plots = 0
    while(os.path.isfile(save_path)):
        save_path = '../Poster_plots/junk/test_plot({})'.format(repeated_plots)
        repeated_plots+=1
    
    print('Saving plot in {}'.format(save_path))
    plt.savefig(save_path)
    fig.close()

def plot_eject():
    name = '281122_MJ1.0_MJ0.1/'
    dir = all_data_dir + name + 'out/'
    data = pluto.Pluto(dir)
    
    (
        time_dat_2,
        a_dat_2,
        e_dat_2,
        period_dat_2,
        mass_dat_2,
    ) = Data_parser_helper.getNbodyInformation_dat(name, 2)
    (
        time_dat_3,
        a_dat_3,
        e_dat_3,
        period_dat_3,
        mass_dat_3,
    ) = Data_parser_helper.getNbodyInformation_dat(name, 3)
    i_min_2, i_max_2 = tools.time_split(time_dat_2, t_min, t_max)
    i_min_3, i_max_3 = tools.time_split(time_dat_3, t_min, t_max)
    gap_time = []
    gap_e = []
    gap_a = []
    for j in range(n_min, n_max+1):
        gap_time.append(j*nts)
        var_data = data.primitive_variable(var, j)[0,:,:]
        R, Phi = np.meshgrid(data.grid['faces'][0], data.grid['faces'][1])
        (
            a_j,
            e_j,
            _, _, _, _, _
        ) = gap.ellipse_find(R, Phi, var_data)
        gap_e.append(abs(e_j))
        gap_a.append(abs(a_j))
    fig, ax = plt.subplots(1,1)
    plotp.poster()
    #plt.grid()
    l1 = ax.plot(time_dat_2[i_min_2:i_max_2], a_dat_2[i_min_2:i_max_2], 'b')
    l2 = ax.plot(time_dat_3[i_min_3:i_max_3], a_dat_3[i_min_3:i_max_3], 'r')
    l3 = ax.plot([t_min,t_max], [a_crit,a_crit], 'k--')
    l4 = ax.plot([t_min,t_max], [kepb,kepb], 'g--')
    l5 = ax.plot(gap_time, gap_a, 'm:')
    ax.set_ylabel(r'Semi-Major Axis ($a_{bin}$)', size = 14)
    ax.set_xlabel(r'Time ($T_{bin}$)', size = 14)
    plt.legend(['Jupiter analogue', 'Neptune Analogue','instability limit', 'Kep47b a', 'gap radius'])
    save_path = '../Poster_plots/Finished_plot_4.png'
    repeated_plots = 0
    while (os.path.isfile(save_path)):
        save_path = '../Poster_plots/Finished_plot_4({}).png'.format(repeated_plots)
        repeated_plots += 1
    
    print('saving plot in {}'.format(save_path))
    fig.savefig(save_path, bbox_inches = 'tight',dpi = 320)
    plt.close(fig)

if __name__ == '__main__':
    plotters = [plot_141122_vs,plot_4_N_dat,plot_4_2ds, plot_4_2ds_kms, plot_eject]
    func_index = UI_helper.selectFunctionsToRun(plotters)
    eval('{}()'.format(plotters[func_index].__name__))
