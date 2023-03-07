import pluto
import tools
import scipy.fft as fft
import numpy as np
import matplotlib.pyplot as plt
import os 

import UI_helper
import Navigation_helper
import Data_parser_helper
from Global_variables import *
import plot_params
import Gap_Finder as gap

def trianglewave(n,av):

    if n % 2 == 0:
        even = True
    elif n % 2 == 1:
        even = False
    else: raise(TypeError)

    if even:
        t = np.arange(-n/2+1, n/2+1)
    elif not even:
        t = np.arange(-(n-1)/2, (n+1)/2)
    
    amp = []

    for i in t:
        if abs(i) >= av:
            amp.append(0.0)
        elif abs(i) <= av:
            amp.append(1.0 - abs(i/av))
    
    return t, np.array(amp)

def fourier_fixer(arr):

    arr_len = int(len(arr))

    if arr_len % 2 == 0:
        even = True
    elif arr_len % 2 == 1:
        even = False
    else: raise(TypeError)

    arr_fix = []

    if even:
        for i in range(0,arr_len-1):
            arr_fix.append(arr[arr_len-1-i])
        for j in range(0,arr_len-1):
            arr_fix.append(arr[j+1])
    elif not even:
        for i in range(0,arr_len-1):
            arr_fix.append(arr[arr_len-1-i])
        for j in arr:
            arr_fix.append(j)
    else: raise(TypeError)

    return arr_fix

def migration_with_fourier() -> None:
    data_name = UI_helper.selectDataToPlot()
    directories = Navigation_helper.Directories(data_name)
    object = UI_helper.selectObjectToPlot(directories.out_dir)
    av = UI_helper.select_averaging_length()

    (
        time,
        a,
        e,
        anomoly,
        mass,
    ) = Data_parser_helper.getNbodyInformation_dat(data_name, object.id)

    delta_a_list = np.array([a[i+1] - a[i] for i in range(len(a)-1)])
    a_delta_t_list = np.array([a[i]*(time[i+1] - time[i]) for i in range(len(a)-1)])
    list_len = int(len(delta_a_list))
    a_dot_over_a_list = delta_a_list/a_delta_t_list
    time_list = np.array(time[0:list_len])

    a_trans_rough = fft.rfft(a_dot_over_a_list)

    _, tri = trianglewave(list_len, av)

    tri_trans_rough = fft.rfft(tri)

    short_len = int(len(a_trans_rough))

    conv = []
    for i in range(0,short_len):
        conv.append(tri_trans_rough[i] * a_trans_rough[i])
    conv = np.array(conv)

    filtered_data = fft.irfft(conv)
    filtered_data_real = np.real(filtered_data)
    filtered_data_real = fft.ifftshift(filtered_data_real)

    filtered_len = int(len(filtered_data_real))

    fig = plt.figure()
    plt.plot(time_list[0:filtered_len-1], filtered_data_real[0:filtered_len-1])
    plt.xlabel('Time (binary orbits)')
    plt.ylabel('a_dot/a (a/binary orbit time)')
    plt.grid()
    fig.tight_layout()

    save_path = '{}obj{}_migration_fourier{}.png'.format(directories.plots_dir, object.id, av)
    repeated_plots = 0
    while (os.path.isfile(save_path)):
        save_path = '{}obj{}_migration_fourier{}({}).png'.format(directories.plots_dir, object.id, av, repeated_plots)
        repeated_plots += 1
    
    print('Saving plot in {0}'.format(save_path))
    fig.savefig(save_path)
    plt.close(fig)

if __name__ == '__main__':
    plotters = [migration_with_fourier]
    func_index = UI_helper.selectFunctionsToRun(plotters)
    eval('{}()'.format(plotters[func_index].__name__))