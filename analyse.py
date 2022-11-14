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
TwoD_plotter = __import__('2d_sigma')
OneD_profile_plotter = __import__('1d_profile')
from Global_variables import *

if __name__ == '__main__':
    plotters = [
        TwoD_plotter.plot_one, 
        TwoD_plotter.plot_many, 
        TwoD_plotter.animate,
        OneD_profile_plotter.plot_one, 
        OneD_profile_plotter.plot_many, 
        OneD_profile_plotter.animate,
        ]
    func_index = UI_helper.selectFunctionsToRun(plotters)
    eval('{}.{}()'.format(plotters[func_index].__module__, plotters[func_index].__name__))