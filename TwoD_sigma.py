#!/usr/bin/env python3
'''
- data.time gives the info about the times - good for iterating over
- data.last gives the last datafile number as int
'''

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

cart = False # 'cart grid'
logsc = False # log scale
nbody = True # nbody integrater used
nts   = 100 # output time step
var   = 'rho' # 'rho, vx1, vx2, prs'

def plot_many() -> None:
    data_name = UI_helper.selectDataToPlot()
    out_dir = all_data_dir + data_name + '/out'
    many_data_files_to_plot = UI_helper.selectManyDataFilesToPlot(out_dir)
    for data_file_to_plot in many_data_files_to_plot:
        plot(data_name, data_file_to_plot)

def plot_one() -> None:
    data_name = UI_helper.selectDataToPlot()
    out_dir = all_data_dir + data_name + '/out'
    data_file_to_plot = UI_helper.selectDataFileToPlot(out_dir)
    plot(data_name, data_file_to_plot)


def plot(data_name: str, data_file_to_plot: int) -> None:

    n = data_file_to_plot
    fig, ax = plt.subplots(1, 1, subplot_kw=dict(aspect='equal',
                                                xlim=[-size,size],
                                                ylim=[-size,size]
                                                ))
    data_parent_dir = all_data_dir + data_name
    out_dir = all_data_dir + data_name + '/out'
    plots_dir = data_parent_dir + '/Plots/'

    data = pluto.Pluto(out_dir)
    sigma = data.primitive_variable(var, n)[0,:,:] #* data.units['density']
    print(np.argwhere(np.isnan(sigma)))
    #sigma[np.isnan(sigma)] = 1.0
    '''
    data.grid is a dictionary of items
    - total: int. the number of datapoints
    - cells: 3 elem list of ints. list of the dimensions of the coord system.
    - faces: 3 elem list of lists. list of positions of the faces of the cells
    - centers: 3 elem list of lists. list of the positions of the centers of the cells
    - dx: 2d list. list of the size of a cell?
    - dV: 2d list. list of the volume of a cell? But all the vals are negative??
    '''

    r, phi = data.grid['centers'].X1, data.grid['centers'].X2

    #ax.set_ylabel(r'$y\;[au]$',fontsize=fs)
    #ax.set_xlabel(r'$x\;[au]$',fontsize=fs)
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
        sig_cb = np.where(r2>1.5, sigma.flatten(), 0.0)
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
        x_0, y_0, phi_0, a, b, num_e = tools.calc_gap_ellipse(sigma, r, phi, fraction=0.1)
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
    if logsc:
        if var=='rho':
            plot = ax.pcolor(X*a_bin, Y*a_bin, np.log10(sigma/np.max(sigma)), cmap='gist_heat')#, vmin=-2)
        elif var=='vx2':
            plot = ax.pcolor(X*a_bin, Y*a_bin, np.log10(sigma*np.sqrt(R[:-1,:-1])/corr), cmap='bwr', vmin=-0.05, vmax=0.05)
        else:
            plot = ax.pcolor(X*a_bin, Y*a_bin, sigma, cmap='bwr', vmin=-0.05, vmax=0.05)
    else:
        plot = ax.pcolor(X*a_bin, Y*a_bin, sigma/np.max(sigma), cmap='viridis')#, vmin=0.0, vmax=0.1)#, vmin=0.01, vmax=2000)#,
    #plot = ax.pcolor(R, phi, sigma, cmap='magma', vmin=1.0, vmax=100.)#, norm=colors.LogNorm(vmin=sigma.min(), vmax=sigma.max()))

    # Colorbar
    #cb = plt.colorbar(plot, orientation='vertical')
    #cb.set_label(r'$log\Sigma\;\left[\mathrm{g}\,/\mathrm{cm}^{-2}\right]$',fontsize=fs)

    # Centre of mass
    ax.plot(0,0, '+k', ms=3)

    ax.text(0.95, 0.95, r't={} yr'.format(n*nts),
            #0.95, 0.95, r't={}$T_\mathrm{{b}}$'.format(n*nts),
            color='c',
            size='x-large',
            #weight='bold',
            ha='right',
            va='top',transform=ax.transAxes)

    # Plot ellipse
    ax.plot(x_ell*a_bin, y_ell*a_bin, '--w', lw=1.5)

    ax.text(0.05, 0.15,
            r'$e_\mathrm{{gap}} = {0:.2f}$'.format(num_e)+'\n'
            #+r'$a_\mathrm{{gap}} = {0:.2f}\,au$'.format(a*a_bin),
            +r'$a_\mathrm{{gap}} = {0:.2f}\,a_\mathrm{{b}}$'.format(a),
            color='w',
            size='x-large',
            #weight='bold',
            ha='left',
            va='top',
            transform=ax.transAxes)

    if nbody:
        # plot the planet evolution here
        txn, pid, x, y = np.loadtxt('{}/nbody.out'.format(out_dir), usecols=(0,1,3,4), unpack=True)
        tim, pid2, aa, ee, varpi = np.loadtxt('{}/nbody_orbital_elements.out'.format(out_dir), usecols=(0,1,2,3,6), unpack=True)
        # Plot Nbody particles
        N=len(np.unique(pid))
        ti = tim[pid2==1]
        tx = txn[pid ==0]
        print('number of bodies = ', N)
        for nn in range(N):
            print(nn)
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
                
                ax.text(0.5, 0.15*nn,
                    r'$e_\mathrm{{p}} = {0:.3e}$'.format(float(e[n]))+'\n'
                    #+r'$a_\mathrm{{gap}} = {0:.2f}\,au$'.format(a*a_bin),
                    +r'$a_\mathrm{{p}} = {0:.3e}\,a_\mathrm{{b}}$'.format(float(a[n])),
                    color='c',
                    size='x-large',
                    #weight='bold',
                    ha='left',
                    va='top',
                    transform=ax.transAxes)
                
    fig.tight_layout() 

    save_path = '{}{}_2d_sigma_{}.png'.format(plots_dir, data_name, n)
    repeated_plots = 1
    while(os.path.isfile(save_path)):
        save_path = '{}{}_2d_sigma_{}({}).png'.format(plots_dir, data_name, n, repeated_plots)
        repeated_plots += 1
    # out_name = '2d_sigma_{}.png'.format(data_name)
    # outfig = out_path+out_name
    print('Saving plot in {0}'.format(save_path))
    fig.savefig(save_path)
    plt.close(fig)

def animate():
    data_name = UI_helper.selectDataToPlot()
    out_dir = all_data_dir + data_name + '/out'
    n_max = Navigation_helper.findMaxFileName(out_dir)

    data_parent_dir = all_data_dir + data_name
    out_dir = all_data_dir + data_name + '/out' 
    plots_dir = data_parent_dir + '/Plots/'

    fig, ax = plt.subplots(1, 1, subplot_kw=dict(aspect='equal',
                                                xlim=[-size,size],
                                                ylim=[-size,size]
                                                ))
    camera = Camera(fig)
    for n in range(n_max+1):
        data = pluto.Pluto(out_dir)
        sigma = data.primitive_variable(var, n)[0,:,:] #* data.units['density']
        # print(np.argwhere(np.isnan(sigma)))
        #sigma[np.isnan(sigma)] = 1.0
        '''
        data.grid is a dictionary of items
        - total: int. the number of datapoints
        - cells: 3 elem list of ints. list of the dimensions of the coord system.
        - faces: 3 elem list of lists. list of positions of the faces of the cells
        - centers: 3 elem list of lists. list of the positions of the centers of the cells
        - dx: 2d list. list of the size of a cell?
        - dV: 2d list. list of the volume of a cell? But all the vals are negative??
        '''

        r, phi = data.grid['centers'].X1, data.grid['centers'].X2

        #ax.set_ylabel(r'$y\;[au]$',fontsize=fs)
        #ax.set_xlabel(r'$x\;[au]$',fontsize=fs)
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
            sig_cb = np.where(r2>1.5, sigma.flatten(), 0.0)
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
            x_0, y_0, phi_0, a, b, num_e = tools.calc_gap_ellipse(sigma, r, phi, fraction=0.1)
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
        if logsc:
            if var=='rho':
                plot = ax.pcolor(X*a_bin, Y*a_bin, np.log10(sigma/np.max(sigma)), cmap='gist_heat')#, vmin=-2)
            elif var=='vx2':
                plot = ax.pcolor(X*a_bin, Y*a_bin, np.log10(sigma*np.sqrt(R[:-1,:-1])/corr), cmap='bwr', vmin=-0.05, vmax=0.05)
            else:
                plot = ax.pcolor(X*a_bin, Y*a_bin, sigma, cmap='bwr', vmin=-0.05, vmax=0.05)
        else:
            plot = ax.pcolor(X*a_bin, Y*a_bin, sigma/np.max(sigma), cmap='viridis')#, vmin=0.0, vmax=0.1)#, vmin=0.01, vmax=2000)#,
        #plot = ax.pcolor(R, phi, sigma, cmap='magma', vmin=1.0, vmax=100.)#, norm=colors.LogNorm(vmin=sigma.min(), vmax=sigma.max()))

        # Colorbar
        #cb = plt.colorbar(plot, orientation='vertical')
        #cb.set_label(r'$log\Sigma\;\left[\mathrm{g}\,/\mathrm{cm}^{-2}\right]$',fontsize=fs)

        # Centre of mass
        ax.plot(0,0, '+k', ms=3)

        ax.text(0.95, 0.95, r't={} yr'.format(n*nts),
                #0.95, 0.95, r't={}$T_\mathrm{{b}}$'.format(n*nts),
                color='c',
                size='x-large',
                #weight='bold',
                ha='right',
                va='top',transform=ax.transAxes)

        # Plot ellipse
        ax.plot(x_ell*a_bin, y_ell*a_bin, '--w', lw=1.5)

        ax.text(0.05, 0.15,
                r'$e_\mathrm{{gap}} = {0:.2f}$'.format(num_e)+'\n'
                #+r'$a_\mathrm{{gap}} = {0:.2f}\,au$'.format(a*a_bin),
                +r'$a_\mathrm{{gap}} = {0:.2f}\,a_\mathrm{{b}}$'.format(a),
                color='w',
                size='x-large',
                #weight='bold',
                ha='left',
                va='top',
                transform=ax.transAxes)

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
                    
                    ax.text(0.5, 0.15*nn,
                        r'$e_\mathrm{{p}} = {0:.3e}$'.format(float(e[n]))+'\n'
                        #+r'$a_\mathrm{{gap}} = {0:.2f}\,au$'.format(a*a_bin),
                        +r'$a_\mathrm{{p}} = {0:.3e}\,a_\mathrm{{b}}$'.format(float(a[n])),
                        color='c',
                        size='x-large',
                        #weight='bold',
                        ha='left',
                        va='top',
                        transform=ax.transAxes)
                    
        fig.tight_layout() 
        # print(' finished plotting {0}'.format(data_name))
        camera.snap()

    save_path = '{}{}_2d_sigma_ANIMATION.gif'.format(plots_dir, data_name)
    repeated_plots = 0
    while(os.path.isfile(save_path)):
        save_path = '{}{}_2d_sigma_ANIMATION({}).png'.format(plots_dir, data_name, repeated_plots)
        repeated_plots += 1

    animation = camera.animate()
    animation.save(save_path)

if __name__ == '__main__':
    plotters = [plot_one, plot_many, animate]
    func_index = UI_helper.selectFunctionsToRun(plotters)
    eval('{}()'.format(plotters[func_index].__name__))