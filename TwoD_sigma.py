# TODO: add pericentre to the legend plot, phi_0

import pluto
import tools
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib as mpl
import matplotlib.colors as colors
import argparse as argp
import os
import Gap_Finder as gap

import UI_helper
import Navigation_helper
from celluloid import Camera
from Global_variables import *
import plot_params
import plotter_helper
from tools import Unit_conv

def plot_many() -> None:
    plot_params.TwoD_sigma_viva()
    data_name = UI_helper.selectDataToPlot()
    out_dir = all_data_dir + data_name + '/out'
    many_data_files_to_plot = UI_helper.selectManyDataFilesToPlot(out_dir)
    for data_file_to_plot in many_data_files_to_plot:
        plot(data_name, data_file_to_plot)

def plot_one() -> None:
    plot_params.TwoD_sigma_viva()
    data_name = UI_helper.selectDataToPlot()
    out_dir = all_data_dir + data_name + '/out'
    data_file_to_plot = UI_helper.selectDataFileToPlot(out_dir)
    plot(data_name, data_file_to_plot)

def plot_one_velocities() -> None:
    data_name = UI_helper.selectDataToPlot()
    out_dir = all_data_dir + data_name + '/out'
    data_file_to_plot = UI_helper.selectDataFileToPlot(out_dir)
    plot_velocities(data_name, data_file_to_plot)

def plot_velocities(data_name: str, data_file_to_plot: int):
    n = data_file_to_plot
    fig, ax = plt.subplots(1, 2, subplot_kw=dict(aspect='equal',
                                                xlim=[-size,size],
                                                ylim=[-size,size]
                                                ))
    directories = Navigation_helper.Directories(data_name)
    data = pluto.Pluto(directories.out_dir)

    global logsc
    global var

    logsc_temp = logsc
    var_temp = var

    logsc = False
    var   = "vx1"
    plot_the_data(n, data, directories.out_dir, ax[0])
    var  = "vx2"
    plot_the_data(n, data, directories.out_dir, ax[1])
                
    fig.tight_layout() 

    log_str = ''
    fname = '{}2d_{}_{}'.format(log_str, var, n)
    save_path = plotter_helper.define_save_plot(directories.plots_dir, fname)
    fig.savefig(save_path)
    plt.close(fig)

    logsc = logsc_temp
    var = var_temp

def plot(data_name: str, data_file_to_plot: int) -> None:
    n = data_file_to_plot
    fig, ax = plt.subplots(1, 1, subplot_kw=dict(aspect='equal',
                                                xlim=[-size,size],
                                                ylim=[-size,size]
                                                ))

    directories = Navigation_helper.Directories(data_name)
    data = pluto.Pluto(directories.out_dir)
    plot_the_data(n, directories.out_dir, data, ax)
                
    fig.tight_layout() 

    if logsc: log_str = 'LOG_'
    else: log_str = ''
    fname = '{}2d_{}_{}'.format(log_str, var, n)
    save_path = plotter_helper.define_save_plot(directories.plots_dir, fname)
    fig.savefig(save_path)
    plt.close(fig)

def animate() -> None:
    data_name = UI_helper.selectDataToPlot()
    directories = Navigation_helper.Directories(data_name)
    data = pluto.Pluto(directories.out_dir)
    n_min, n_max = UI_helper.selectAnimateRange(directories.out_dir)

    fig, ax = plt.subplots(1, 1, subplot_kw=dict(aspect='equal',
                                                xlim=[-size,size],
                                                ylim=[-size,size]
                                                ))

    camera = Camera(fig)

    # animate = plt_anim.FuncAnimation(fig, partial(plot_the_data, out_dir=out_dir, ax=ax), n_max)

    for n in range(n_min, n_max+1):
        plot_the_data(n, directories.out_dir, data, ax, show_colorbars=False)
        fig.tight_layout() 
        camera.snap()
        # cb.remove()

    if logsc: log_str = 'LOG_'
    else: log_str = ''
    fname = '{}2d_{}_ANIMATION_{}-{}'.format(log_str, var, n_min, n_max)
    save_path = plotter_helper.define_save_plot(directories.plots_dir, fname, 'gif')

    # animate.save(save_path)

    animation = camera.animate(interval=100)
    animation.save(save_path)
    plt.close(fig)

def plot_the_data(
        n: int, 
        out_dir: str, 
        data: pluto.Pluto, 
        ax: plt.Axes, 
        show_colorbars: bool = True, 
        show_meta_data = True,
        show_instability_zone = False,
        show_Kepler_47_planets = False,
        show_contours = False,
        ):
    var_data = data.primitive_variable(var, n)[0,:,:] #* data.units['density']
    a_bin_in_distance_unit = Unit_conv.distance(a_bin)

    r, phi = data.grid['centers'].X1, data.grid['centers'].X2

    # Circle around cylinder, the large central cell
    r_cylinder = data.grid['faces'][0][0]
    cylinder = patches.Arc((0,0), width=2*Unit_conv.distance(r_cylinder),
                                height=2*Unit_conv.distance(r_cylinder))
    ax.add_patch(cylinder)

    # These are both 2D meshes, but R[φ][r] is the same for all θ, Phi[φ][r] is the same for all r
    R, Phi = np.meshgrid(data.grid['faces'][0], data.grid['faces'][1])

    # Calculation of gap ellipse
    if cart:
        X=R[:-1,:-1]
        Y=Phi[:-1,:-1]
        frac=0.1
        
        r2  = (X*X+Y*Y).flatten()
        phi = np.arctan2(Y,X).flatten()
        sig_cb = np.where(r2>1.5, var_data.flatten(), 0.0)
        apo    = r2[np.argmax(sig_cb)]**0.5
        phi_0  = phi[np.argmax(sig_cb)]
        n_phi_apo = np.argwhere(np.abs(phi-phi_0)<1e-3)
        n_phi_per = np.argwhere(np.abs(phi+phi_0)<1e-3)
        n_per = np.argmax(sig_cb[n_phi_per])
        peri  = r2[n_phi_per[n_per]]**0.5
        
        n_phi_al = np.argwhere((np.abs(phi-phi_0)<1e-3) & (r2<=apo**2))
        n_phi_pl = np.argwhere((np.abs(phi+phi_0)<1e-3) & (r2<=peri**2))
        n_afrac = np.argmin(np.abs(np.max(sig_cb[n_phi_apo])*frac-sig_cb[n_phi_al]))
        n_pfrac = np.argmin(np.abs(np.max(sig_cb[n_phi_per])*frac-sig_cb[n_phi_pl]))

        #print(phi[n_phi_pl[n_pfrac]], phi[n_phi_per[n_per]], peri, r2[n_phi_pl[n_pfrac]]**0.5, phi_0)

        peri_frac   = r2[n_phi_pl[n_pfrac]]**0.5
        apo_frac    = r2[n_phi_al[n_afrac]]**0.5

        a = ((apo_frac + peri_frac)/2.0)[0]
        ea = apo_frac - a
        num_e = (apo_frac - a)/a
        num_e = float(num_e)
        b = (a**2-ea**2)**0.5
        x_0 = ea*np.cos(phi_0)
        y_0 = ea*np.sin(phi_0)
        
        x_ell, y_ell = tools.get_ellipse_points(x_0, y_0, phi_0, a, b)
    else:
        X = np.multiply(R, np.cos(Phi))
        Y = np.multiply(R, np.sin(Phi))
        # a, num_e, phi_0, x_0, y_0, x_ell, y_ell = gap.ellipse_find(R, Phi, var_data)
        x_0, y_0, phi_0, a, b, num_e = tools.calc_gap_ellipse(var_data, r, phi, fraction=0.1)
        x_ell, y_ell = tools.get_ellipse_points(x_0, y_0, phi_0, a, b)
    '''
    vr = data.primitive_variable('vx1', n)[0,:,:]
    vphi = data.primitive_variable('vx2', n)[0,:,:]
    vx = vr*np.cos(Phi[:-1,:-1])-vphi*np.sin(Phi[:-1,:-1]) 
    vy = vr*np.sin(Phi[:-1,:-1])+vphi*np.cos(Phi[:-1,:-1]) 
    ecc = (1+2*((vr**2+vphi**2)/2.-1/R[:-1,:-1])*(Y[:-1,:-1]*vx-X[:-1,:-1]*vy)**2)**0.5
    '''
    corr = np.sqrt(1+(-1.0-0.5)*0.05*0.05)
    # 2D-Plot
    if var=='vx1':
        # TODO: change so the radial velocity isn't scaled to kepler velocity
        if logsc: print('\nWarning: logsc==True, but vx can be -ve. \nplotting not log\n')
        plot = ax.pcolor(
            X*a_bin_in_distance_unit, 
            Y*a_bin_in_distance_unit, 
            var_data, 
            cmap='bwr'
        ) #, vmin=-0.05, vmax=0.05)
        colourbar_label = var
    elif var=='vx2':
        if logsc: print('\nWarning: logsc==True, but vx can be -ve. \nplotting not log\n')
        # plot = ax.pcolor(X*a_bin, Y*a_bin, (var_data*np.sqrt(R[:-1,:-1])/corr), cmap='bwr') #, vmin=-0.05, vmax=0.05)
        plot = ax.pcolor(
            X*a_bin_in_distance_unit, 
            Y*a_bin_in_distance_unit, 
            np.log(var_data*np.sqrt(R[:-1,:-1])/corr), 
            cmap='bwr', 
            vmin=-0.1, 
            vmax=0.1
        )
        colourbar_label = var
    elif var=='rho':
        if logsc:
            plot = ax.pcolor(
                X*a_bin_in_distance_unit, 
                Y*a_bin_in_distance_unit, 
                np.log10(Unit_conv.surface_density(var_data, 'grams', 'cm')), # /np.max(var_data) 
                cmap='gist_heat',
            ) #, vmin=-3)#, vmin=-2)
            colourbar_label = r'$log\Sigma$' + ' [' + Unit_conv.surface_density_label('grams', 'cm') + ']'
            # fr'$log\Sigma\;\left[{{{Unit_conv.surface_density_label()}}}\right]$'
        else:
            plot = ax.pcolor(
                X*a_bin_in_distance_unit, 
                Y*a_bin_in_distance_unit, 
                Unit_conv.surface_density(var_data, 'grams', 'cm'), 
                cmap='viridis',
            )#, vmin=0.0, vmax=0.1)#, vmin=0.01, vmax=2000)#,
            colourbar_label = r'$\Sigma$' + ' [' + Unit_conv.surface_density_label('grams', 'cm') + ']' #fr'$log\Sigma\;\left[{{{Unit_conv.surface_density_label()}}}\right]$'
    else:
        raise Exception('invalid var chosen')

    if show_contours:
        max_density = np.max(Unit_conv.surface_density(var_data, 'grams', 'cm'))
        min_density = np.min(Unit_conv.surface_density(var_data, 'grams', 'cm'))
        levels = min_density + ( (np.linspace(0.0, 1.0, num=10)**2) * max_density )
        contours_plot = plt.contour(
            [X[i][:-1] for i in range(len(X)-1)]*a_bin_in_distance_unit, 
            [Y[i][:-1] for i in range(len(X)-1)]*a_bin_in_distance_unit, 
            Unit_conv.surface_density(var_data, 'grams', 'cm'),
            colors = 'blue',
            levels = levels,
            linewidths=0.5,
        )

    # if logsc:
    #     if var=='rho':
    #         plot = ax.pcolor(X*a_bin, Y*a_bin, np.log10(var_data/np.max(var_data)), cmap='gist_heat')#, vmin=-2)
    #         colourbar_label = r'$log\Sigma\;\left[\mathrm{g}\,/\mathrm{cm}^{-2}\right]$'
    #     elif var=='vx1' or var=='vx2':
    #         print('\nwarning: plotting velocity, which may be -ve, with log scale\n')
    #         plot = ax.pcolor(X*a_bin, Y*a_bin, np.log10(var_data*np.sqrt(R[:-1,:-1])/corr), cmap='bwr', vmin=-0.05, vmax=0.05)
    #         colourbar_label = 'log(velocity)'
    #     else:
    #         raise Exception('invalid var chosen')
    #         # plot = ax.pcolor(X*a_bin, Y*a_bin, var_data, cmap='bwr', vmin=-0.05, vmax=0.05)
    #         # colourbar_label = 'unknown?'
    # else:
    #     if var=='rho':
    #         plot = ax.pcolor(X*a_bin, Y*a_bin, var_data/np.max(var_data), cmap='viridis')#, vmin=0.0, vmax=0.1)#, vmin=0.01, vmax=2000)#,
    #         colourbar_label = r'$log\Sigma\;\left[\mathrm{g}\,/\mathrm{cm}^{-2}\right]$'
    #     elif var=='vx1' or var=='vx2':
    #         plot = ax.pcolor(X*a_bin, Y*a_bin, var_data*np.sqrt(R[:-1,:-1])/corr, cmap='bwr', vmin=-0.05, vmax=0.05)
    #         colourbar_label = 'velocity'
    #     else:
    #         raise Exception('invalid var chosen')

    #plot = ax.pcolor(R, phi, sigma, cmap='magma', vmin=1.0, vmax=100.)#, norm=colors.LogNorm(vmin=sigma.min(), vmax=sigma.max()))

    ax.set_ylabel('y [' + Unit_conv.distance_label() + ']',fontsize=fs)
    ax.set_xlabel('x [' + Unit_conv.distance_label() + ']',fontsize=fs)

    if show_meta_data:
        time_text = ax.text(
            0.95, 
            0.95, 
            r'time={} {}'.format(Unit_conv.time(n*nts), Unit_conv.time_label()),
            #0.95, 0.95, r't={}$T_\mathrm{{b}}$'.format(n*nts),
            color='w',
            size='x-large',
            #weight='bold',
            ha='right',
            va='top',
            transform=ax.transAxes,
            bbox = (dict(facecolor='red', alpha=0.5, edgecolor='red')),
        )
        # time_text.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='black'))

        binary_text = ax.text(
            0.05, 
            0.15,
            r'$e_\mathrm{{gap}} = {0:.2f}$'.format(num_e)+'\n'
            #+r'$a_\mathrm{{gap}} = {0:.2f}\,au$'.format(a*a_bin),
            +r'$a_\mathrm{{gap}} = {0:.2f}\,$'.format(Unit_conv.distance(a))
            +Unit_conv.distance_label(),
            color='w',
            size='x-large',
            #weight='bold',
            ha='left',
            va='top',
            transform=ax.transAxes,
            bbox = dict(facecolor='red', alpha=0.5, edgecolor='red'),
        )
        # Plot ellipse
        ax.plot(
            x_ell*a_bin_in_distance_unit, 
            y_ell*a_bin_in_distance_unit, 
            '--w', 
            lw=1.5
        )
    # Colorbar
    if show_colorbars:
        cb = plt.colorbar(plot, orientation='vertical')
        cb.set_label(colourbar_label, fontsize=fs)
        if show_contours: 
            cb.add_lines(contours_plot)

    
    # Centre of mass
    ax.plot(0,0, '+k', ms=3)

    # Instability Zone
    if show_instability_zone:
        plotter_helper.plot_instability_zone_for_twoD_sigma(ax)
    
    # Kepler-47 Planets
    if show_Kepler_47_planets:
        plotter_helper.plot_Kepler_47_planets_for_twoD_sigma(ax)

    nbody_text_list = []
    if nbody:
        # plot the planet evolution here
        txn, pid, x, y = np.loadtxt('{}/nbody.out'.format(out_dir), usecols=(0,1,3,4), unpack=True)
        tim, pid2, aa, ee, varpi = np.loadtxt('{}/nbody_orbital_elements.out'.format(out_dir), usecols=(0,1,2,3,6), unpack=True)
        # Plot Nbody particles
        N=len(np.unique(pid))
        ti = tim[pid2==1]
        tx = txn[pid ==0]
        for nn in range(N):
            xp, yp = x[pid==nn], y[pid==nn]
            ax.plot(xp[tx==n]*a_bin_in_distance_unit, yp[tx==n]*a_bin_in_distance_unit, '.b', ms=6.0)

            if show_meta_data and nn>1: 
                # if nn>0:
                #print(e[pid2==1], pid2[pid2==nn])
                a = aa[pid2==nn]
                e = ee[pid2==nn]
                #varpi = varpi[pid2==2][n]
                phi = np.arange(0.0, 2.0*np.pi, 0.01)
                r = a[n]*(1.0-e[n]**2)/(1.0-e[n]*np.cos(phi-phi_0))
                x_ell = r*np.cos(phi)
                y_ell = r*np.sin(phi)
                ax.plot(
                    Unit_conv.distance(x_ell),
                    Unit_conv.distance(y_ell),
                    'w:'
                )
                
                e_id = r'$e_{{{}}}$'.format(objects_in_kep47[nn].shorthand_name)
                e_num = r'${0:.3f}$'.format(e[n])
                a_id = r'$a_{{{}}}$'.format(objects_in_kep47[nn].shorthand_name)
                a_num = r'${0:.3f}$'.format(a[n])
            
                nbody_text = ax.text(
                    0.5, 
                    -0.15+0.15*nn,
                    e_id + r'=' + e_num +'\n' + a_id + r'=' + a_num + Unit_conv.distance_label(),
                    color='w',
                    size='x-large',
                    #weight='bold',
                    ha='left',
                    va='top',
                    transform=ax.transAxes,
                    bbox = dict(facecolor='red', alpha=0.5, edgecolor='red')
                )


                nbody_text_list.append(nbody_text)
    if show_meta_data:
        return time_text, binary_text, *nbody_text_list

# def plot_one():
#     n: int
#     fig = plot
#     ...
#     plt.savefig()

# def plot_many():
#     n_list = [n1, n2, ...]
#     for n in n_list:
#         fig = plot
#         ...
#         plt.savefig()

# def animate():
#     n_list = [n1, n2, ...]
#     fig = plot
#     for n in n_list:
#         ...
#         camera.snap(fig)



if __name__ == '__main__':
    plotters = [plot_one, plot_many, plot_one_velocities, animate]
    func_index = UI_helper.selectFunctionsToRun(plotters)
    eval('{}()'.format(plotters[func_index].__name__))