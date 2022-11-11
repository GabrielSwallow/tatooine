#!/usr/bin/env python3
'''
- data.time gives the info about the times - good for iterating over
- data.last gives the last datafile number as int
'''

from distutils.log import error
import pluto
import tools
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib as mpl
import matplotlib.colors as colors
import argparse as argp
import os
from celluloid import Camera

quick_cbd_dir = '/Users/gabe/Documents/Physics/MSci/test_case/quick_cbd/'

datasets = [
    '1_Orbit',
    '10_Orbits',
    '1000_Orbits',
    '1000_Orbits_full'
]

for d_index, filename in enumerate(datasets):
    print(d_index, ' : ', filename)
dataset_index = int(input("please choose a dataset to plot for \n"))
data_name = datasets[dataset_index]

data_dir = quick_cbd_dir + data_name

def findMaxFileName(dataset: str) -> int:
    files = os.listdir(data_dir)
    dataFiles = list(filter(lambda f: f[0:4]=='data', files))
    numbers = [int(d[5:9]) for d in dataFiles]
    maxFileNumber = max(numbers)
    print('file to plot =', max(numbers))
    return maxFileNumber

def findNumBodies(dataset: str) -> int:
    _, pid, _, _ = np.loadtxt( 'nbody.out'.format(data_dir), usecols=(0,1,3,4), unpack=True)
    return len(np.unique(pid))


parser = argp.ArgumentParser(description = '2d surface density plot from file "--number"')
parser.add_argument('--number', default = 0, type = int, help = 'File number')
parser.add_argument('--cart', default = False, type = bool, help = 'cart grid')
parser.add_argument('--log', default = False, type = bool, help = 'log scale')
parser.add_argument('--nbody', default = False, type = bool, help = 'nbody integrater used')
parser.add_argument('--size', default = 7.5, type = float, help = 'size of the plot')
parser.add_argument('--time', default = 100, type = float, help = 'output time step')
parser.add_argument('--var', default = 'rho', type = str, help = 'rho, vx1, vx2, prs')
args = parser.parse_args()

a_bin = 1
size = args.size
Rmax = 70
# File number
num_data_files = findMaxFileName(data_name)
cart = args.cart
logsc = args.log
nbody = True
nts   = args.time
var   = args.var
# Number of nbodies
fs = 14 # font size

C0, C1, C2 = 'k', 'b', 'y'

size = size*a_bin
def plot(datafile: int):
    n = datafile
    fig, ax = plt.subplots(1, 1, subplot_kw=dict(aspect='equal',
                                                xlim=[-size,size],
                                                ylim=[-size,size]
                                                ))

    data = pluto.Pluto(data_dir)
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
        txn, pid, x, y = np.loadtxt('{}/nbody.out'.format(data_dir), usecols=(0,1,3,4), unpack=True)
        tim, pid2, aa, ee, varpi = np.loadtxt('{}/nbody_orbital_elements.out'.format(data_dir), usecols=(0,1,2,3,6), unpack=True)
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

    plots_dir = '/Users/gabe/Documents/Physics/MSci/test_case/Plots/'

    save_path = '{}{}_2d_sigma.png'.format(plots_dir, data_name)
    repeated_plots = 0
    while(os.path.isfile(save_path)):
        save_path = '{}{}_2d_sigma_{}.png'.format(plots_dir, data_name, repeated_plots)
        repeated_plots += 1
    # out_name = '2d_sigma_{}.png'.format(data_name)
    # outfig = out_path+out_name
    print('Saving plot in {0}'.format(save_path))
    fig.savefig(save_path)
    plt.close(fig)

def animate():
    fig, ax = plt.subplots(1, 1, subplot_kw=dict(aspect='equal',
                                                xlim=[-size,size],
                                                ylim=[-size,size]
                                                ))
    camera = Camera(fig)
    for n in range(num_data_files):
        data = pluto.Pluto(data_dir)
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
            txn, pid, x, y = np.loadtxt('{}/nbody.out'.format(data_dir), usecols=(0,1,3,4), unpack=True)
            tim, pid2, aa, ee, varpi = np.loadtxt('{}/nbody_orbital_elements.out'.format(data_dir), usecols=(0,1,2,3,6), unpack=True)
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

        plots_dir = '/Users/gabe/Documents/Physics/MSci/test_case/Plots/'

        save_path = '{}{}_2d_sigma.png'.format(plots_dir, data_name)
        repeated_plots = 0
        while(os.path.isfile(save_path)):
            save_path = '{}{}_2d_sigma_{}.png'.format(plots_dir, data_name, repeated_plots)
            repeated_plots += 1
        # out_name = '2d_sigma_{}.png'.format(data_name)
        # outfig = out_path+out_name
        print(' finished plotting {0}'.format(data_name))
        camera.snap()
    animation = camera.animate()
    animation.save('animation.gif')


if __name__ == '__main__':
    # for datafile in range(num_data_files):
    #     plot(datafile)
    
    # animate()
    pass
    