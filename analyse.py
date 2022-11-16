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

import TwoD_sigma
import OneD_profile
import Nbody_Characteristics

plotters = [
    TwoD_sigma.plot_one, 
    TwoD_sigma.plot_many, 
    TwoD_sigma.animate,
    OneD_profile.plot_one, 
    OneD_profile.plot_many, 
    OneD_profile.animate,
    Nbody_Characteristics.plot_one,
]

if __name__ == '__main__':
    ui_input = UI_helper.selectFunctionsToRun(plotters)
    if type(ui_input) == int:
        eval('{}.{}()'.format(plotters[ui_input].__module__, plotters[ui_input].__name__))
    elif ui_input == 'all':
        for f in plotters:
            eval('{}.{}()'.format(f.__module__, f.__name__))
