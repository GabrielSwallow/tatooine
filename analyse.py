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
import TwoD_sigma
import OneD_profile
from Global_variables import *

if __name__ == '__main__':
    plotters = [
        TwoD_sigma.plot_one, 
        TwoD_sigma.plot_many, 
        TwoD_sigma.animate,
        OneD_profile.plot_one, 
        OneD_profile.plot_many, 
        OneD_profile.animate,
        ]
    func_index = UI_helper.selectFunctionsToRun(plotters)
    eval('{}.{}()'.format(plotters[func_index].__module__, plotters[func_index].__name__))