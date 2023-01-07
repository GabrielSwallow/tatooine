import pluto
import tools
import numpy as np
import matplotlib.pyplot as plt
from celluloid import Camera
import os 

import UI_helper
import Navigation_helper
import Data_parser_helper
from Global_variables import *

def plot_one() -> None:
    data_name = UI_helper.selectDataToPlot()
    plot(data_name)

def plot(data_name: str) -> None:
    directories = Navigation_helper.Directories(data_name)
    
    # obj, obj_des = UI_helper.selectObjectToPlot(directories.out_dir)
    (time, _, _, _, _, _, _, _, accretion) = Data_parser_helper.get_averages_data(data_name) 

    fig = plt.figure()
    plt.plot(time, accretion)
    fig.tight_layout()
    
    save_path = '{}{}_obj2_accretion.png'.format(directories.plots_dir, data_name)
    repeated_plots = 0
    while(os.path.isfile(save_path)):
        save_path = '{}{}_obj2__accretion({}).png'.format(directories.plots_dir, data_name, repeated_plots)
        repeated_plots += 1
    
    print('Saving plot in {0}'.format(save_path))
    fig.savefig(save_path)
    plt.close(fig)
    
if __name__ == '__main__':
    plotters = []
    func_index = UI_helper.selectFunctionsToRun(plotters)
    eval('{}()'.format(plotters[func_index].__name__))