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

def plot_many() -> None:
    plot_params.default()
    data_name = UI_helper.selectDataToPlot()
    out_dir = all_data_dir + data_name + '/out'
    many_data_files_to_plot = UI_helper.selectManyDataFilesToPlot(out_dir)
    for data_file_to_plot in many_data_files_to_plot:
        plot(data_name, data_file_to_plot)

def plot_one() -> None:
    plot_params.default()
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
    save_path = '{}{}2d_{}_{}.png'.format(directories.plots_dir, log_str, var, n)
    repeated_plots = 1
    while(os.path.isfile(save_path)):
        save_path = '{}{}2d_{}_{}({}).png'.format(directories.plots_dir, log_str, var, n, repeated_plots)
        repeated_plots += 1
    print('Saving plot in {0}'.format(save_path))
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
    save_path = '{}{}2d_{}_{}.png'.format(directories.plots_dir, log_str, var, n)
    repeated_plots = 1
    while(os.path.isfile(save_path)):
        save_path = '{}{}2d_{}_{}({}).png'.format(directories.plots_dir, log_str, var, n, repeated_plots)
        repeated_plots += 1
    print('Saving plot in {0}'.format(save_path))
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
        plot_the_data(n, directories.out_dir, data, ax, plot_colorbars=False)
        fig.tight_layout() 
        camera.snap()
        # cb.remove()

    if logsc: log_str = 'LOG_'
    else: log_str = ''
    save_path = '{}{}2d_{}_ANIMATION_{}-{}.gif'.format(directories.plots_dir, log_str, var, n_min, n_max)

    repeated_plots = 0
    while(os.path.isfile(save_path)):
        save_path = '{}{}2d_{}_ANIMATION_{}-{}({}).gif'.format(directories.plots_dir, log_str, var, n_min, n_max, repeated_plots)
        repeated_plots += 1

    # animate.save(save_path)

    animation = camera.animate()
    animation.save(save_path)
    plt.close(fig)

def plot_the_data(n: int, out_dir: str, data: pluto.Pluto, ax: plt.Axes, plot_colorbars: bool = True):
    var_data = data.primitive_variable(var, n)[0,:,:] #* data.units['density']

    r, phi = data.grid['centers'].X1, data.grid['centers'].X2

    ax.set_ylabel(r'$y\;[a_\mathrm{b}]$',fontsize=fs)
    ax.set_xlabel(r'$x\;[a_\mathrm{b}]$',fontsize=fs)

    # Circle around cylinder, the large central cell
    r_cylinder = data.grid['faces'][0][0]
    cylinder = patches.Arc((0,0), width=2*r_cylinder,
                                height=2*r_cylinder)
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
        a, num_e, phi_0, x_0, y_0, x_ell, y_ell = gap.ellipse_find(R, Phi, var_data)
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
        plot = ax.pcolor(X*a_bin, Y*a_bin, var_data, cmap='bwr') #, vmin=-0.05, vmax=0.05)
        colourbar_label = var
    elif var=='vx2':
        if logsc: print('\nWarning: logsc==True, but vx can be -ve. \nplotting not log\n')
        # plot = ax.pcolor(X*a_bin, Y*a_bin, (var_data*np.sqrt(R[:-1,:-1])/corr), cmap='bwr') #, vmin=-0.05, vmax=0.05)
        plot = ax.pcolor(X*a_bin, Y*a_bin, np.log(var_data*np.sqrt(R[:-1,:-1])/corr), cmap='bwr', vmin=-0.1, vmax=0.1)
        colourbar_label = var
    elif var=='rho':
        if logsc:
            plot = ax.pcolor(X*a_bin, Y*a_bin, np.log10(var_data/np.max(var_data)), cmap='gist_heat')#, vmin=-2)
            colourbar_label = r'$log\Sigma\;\left[\mathrm{g}\,/\mathrm{cm}^{-2}\right]$'
        else:
            plot = ax.pcolor(X*a_bin, Y*a_bin, var_data/np.max(var_data), cmap='viridis')#, vmin=0.0, vmax=0.1)#, vmin=0.01, vmax=2000)#,
            colourbar_label = r'$log\Sigma\;\left[\mathrm{g}\,/\mathrm{cm}^{-2}\right]$'
    else:
        raise Exception('invalid var chosen')

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

    # Colorbar
    if plot_colorbars:
        cb = plt.colorbar(plot, orientation='vertical')
        cb.set_label(colourbar_label, fontsize=fs)

    # Centre of mass
    ax.plot(0,0, '+k', ms=3)

    time_text = ax.text(0.95, 0.95, r't={} yr'.format(n*nts),
            #0.95, 0.95, r't={}$T_\mathrm{{b}}$'.format(n*nts),
            color='c',
            size='x-large',
            #weight='bold',
            ha='right',
            va='top',transform=ax.transAxes)

    # Plot ellipse
    ax.plot(x_ell*a_bin, y_ell*a_bin, '--w', lw=1.5)

    binary_text = ax.text(0.05, 0.15,
           r'$e_\mathrm{{gap}} = {0:.2f}$'.format(num_e)+'\n'
           #+r'$a_\mathrm{{gap}} = {0:.2f}\,au$'.format(a*a_bin),
           +r'$a_\mathrm{{gap}} = {0:.2f}\,a_\mathrm{{b}}$'.format(a),
           color='w',
           size='x-large',
           #weight='bold',
           ha='left',
           va='top',
           transform=ax.transAxes)
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
            ax.plot(xp[tx==n]*a_bin, yp[tx==n]*a_bin, '.b', ms=6.0)
            if nn>0:
                #print(e[pid2==1], pid2[pid2==nn])
                a = aa[pid2==nn]
                e = ee[pid2==nn]
                #varpi = varpi[pid2==2][n]
                phi = np.arange(0.0, 2.0*np.pi, 0.01)
                r = a[n]*(1.0-e[n]**2)/(1.0-e[n]*np.cos(phi-phi_0))
                x_ell = r*np.cos(phi)
                y_ell = r*np.sin(phi)
                ax.plot(x_ell,y_ell,'w:')
                
                nbody_text = ax.text(0.5, 0.15*nn,
                    r'$e_\mathrm{{p}} = {0:.3e}$'.format(float(e[n]))+'\n'
                    #+r'$a_\mathrm{{gap}} = {0:.2f}\,au$'.format(a*a_bin),
                    +r'$a_\mathrm{{p}} = {0:.3e}\,a_\mathrm{{b}}$'.format(float(a[n])),
                    color='c',
                    size='x-large',
                    #weight='bold',
                    ha='left',
                    va='top',
                    transform=ax.transAxes)
                nbody_text_list.append(nbody_text)
    
    return time_text, binary_text, *nbody_text_list
                    
        
def plot_many_alt() -> None:
    plot_params.default()
    data_name = UI_helper.selectDataToPlot()
    out_dir = all_data_dir + data_name + '/out'
    many_data_files_to_plot = UI_helper.selectManyDataFilesToPlot(out_dir)
    for data_file_to_plot in many_data_files_to_plot:
        plot_alt(data_name, data_file_to_plot)

def plot_alt(data_name: str, data_file_to_plot: int) -> None:
    n = data_file_to_plot
    fig, ax = plt.subplots(1, 1, subplot_kw=dict(aspect='equal',
                                                xlim=[-size,size],
                                                ylim=[-size,size]
                                                ))

    directories = Navigation_helper.Directories(data_name)
    data = pluto.Pluto(directories.out_dir)
    plot_the_data_alt(n, directories.out_dir, data, ax)
                
    fig.tight_layout() 
    plot_params.poster

    if logsc: log_str = 'LOG_'
    else: log_str = ''
    save_path = '{}{}2d_{}_{}.png'.format(directories.plots_dir, log_str, var, n)
    repeated_plots = 1
    while(os.path.isfile(save_path)):
        save_path = '{}{}2d_{}_{}({}).png'.format(directories.plots_dir, log_str, var, n, repeated_plots)
        repeated_plots += 1
    print('Saving plot in {0}'.format(save_path))
    fig.savefig(save_path, bbox_inches = 'tight', dpi = 160)
    plt.close(fig)

def plot_the_data_alt(n: int, out_dir: str, data: pluto.Pluto, ax: plt.Axes, plot_colorbars: bool = True):
    var_data = data.primitive_variable(var, n)[0,:,:] #* data.units['density']

    r, phi = data.grid['centers'].X1, data.grid['centers'].X2

    ax.set_ylabel(r'$y\;[a_\mathrm{b}]$',fontsize=20)
    ax.set_xlabel(r'$x\;[a_\mathrm{b}]$',fontsize=20)

    # Circle around cylinder, the large central cell
    r_cylinder = data.grid['faces'][0][0]
    cylinder = patches.Arc((0,0), width=2*r_cylinder,
                                height=2*r_cylinder)
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
        a, num_e, phi_0, x_0, y_0, x_ell, y_ell = gap.ellipse_find(R, Phi, var_data)
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
        plot = ax.pcolor(X*a_bin, Y*a_bin, var_data, cmap='bwr') #, vmin=-0.05, vmax=0.05)
        colourbar_label = var
    elif var=='vx2':
        if logsc: print('\nWarning: logsc==True, but vx can be -ve. \nplotting not log\n')
        # plot = ax.pcolor(X*a_bin, Y*a_bin, (var_data*np.sqrt(R[:-1,:-1])/corr), cmap='bwr') #, vmin=-0.05, vmax=0.05)
        plot = ax.pcolor(X*a_bin, Y*a_bin, np.log(var_data*np.sqrt(R[:-1,:-1])/corr), cmap='bwr', vmin=-0.1, vmax=0.1)
        colourbar_label = var
    elif var=='rho':
        if logsc:
            plot = ax.pcolor(X*a_bin, Y*a_bin, np.log10(var_data/np.max(var_data)), cmap='gist_heat', vmin=-3)
            colourbar_label = r'$log\Sigma\;\left[\mathrm{g}\,/\mathrm{cm}^{-2}\right]$'
        else:
            plot = ax.pcolor(X*a_bin, Y*a_bin, var_data/np.max(var_data), cmap='viridis')#, vmin=0.0, vmax=0.1)#, vmin=0.01, vmax=2000)#,
            colourbar_label = r'$log\Sigma\;\left[\mathrm{g}\,/\mathrm{cm}^{-2}\right]$'
    else:
        raise Exception('invalid var chosen')

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

    # Colorbar
    if plot_colorbars:
        cb = plt.colorbar(plot, orientation='vertical')
        cb.set_label(colourbar_label, fontsize=fs)

    # Centre of mass
    ax.plot(0,0, '+k', ms=3)

    time_text = ax.text(0.95, 0.98, r't={} yr'.format(n*nts),
            #0.95, 0.95, r't={}$T_\mathrm{{b}}$'.format(n*nts),
            color='b',
            size=20,
            #weight='bold',
            ha='right',
            va='top',transform=ax.transAxes)

    # Plot ellipse
    ax.plot(x_ell*a_bin, y_ell*a_bin, '--w', lw=1.5)

    binary_text = ax.text(0.05, 0.20,
            r'$e_\mathrm{{gap}} = {0:.2f}$'.format(num_e)+'\n'
            #+r'$a_\mathrm{{gap}} = {0:.2f}\,au$'.format(a*a_bin),
            +r'$a_\mathrm{{gap}} = {0:.2f}\,a_\mathrm{{b}}$'.format(a),
            color='w',
            size=20,
            #weight='bold',
            ha='left',
            va='top',
            transform=ax.transAxes)
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
            ax.plot(xp[tx==n]*a_bin, yp[tx==n]*a_bin, '.b', ms=6.0)
            if nn>0:
                #print(e[pid2==1], pid2[pid2==nn])
                a = aa[pid2==nn]
                e = ee[pid2==nn]
                #varpi = varpi[pid2==2][n]
                phi = np.arange(0.0, 2.0*np.pi, 0.01)
                r = a[n]*(1.0-e[n]**2)/(1.0-e[n]*np.cos(phi-phi_0))
                x_ell = r*np.cos(phi)
                y_ell = r*np.sin(phi)
                #ax.plot(x_ell,y_ell,'w:')
                
                # nbody_text = ax.text(0.5, 0.15*nn,
                #     r'$e_\mathrm{{p}} = {0:.3e}$'.format(float(e[n]))+'\n'
                #     #+r'$a_\mathrm{{gap}} = {0:.2f}\,au$'.format(a*a_bin),
                #     +r'$a_\mathrm{{p}} = {0:.3e}\,a_\mathrm{{b}}$'.format(float(a[n])),
                #     color='c',
                #     size='x-large',
                #     #weight='bold',
                #     ha='left',
                #     va='top',
                #     transform=ax.transAxes)
                # nbody_text_list.append(nbody_text)






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
    plotters = [plot_one, plot_many, plot_one_velocities, animate, plot_many_alt]
    func_index = UI_helper.selectFunctionsToRun(plotters)
    eval('{}()'.format(plotters[func_index].__name__))