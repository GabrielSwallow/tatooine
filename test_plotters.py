'''
Simple test to make sure all plotters run, and don't crash
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

import TwoD_sigma
import OneD_profile
import Nbody_Characteristics
import Accretion
import migration
import Torque

from unittest.mock import patch, create_autospec

auto_plotters = [
    # TwoD_sigma.plot_one, 
    TwoD_sigma.plot_many, 
    # TwoD_sigma.animate,
    # OneD_profile.plot_one, 
    OneD_profile.plot_many, 
    # OneD_profile.animate,
    # Nbody_Characteristics.plot_one_using_dat,
    Nbody_Characteristics.plot_one_using_out,
    Accretion.plot_one_disk_accretion,
    Accretion.plot_one_planet_accretion,
    migration.calculate_migration,
    Torque.plot_torque,
]

data_name = UI_helper.selectDataToPlot()
directories = Navigation_helper.Directories(data_name)
n_max = Navigation_helper.findMaxFileNumber(directories.out_dir)

@patch('UI_helper.selectFunctionsToRun')
@patch('UI_helper.selectObjectToPlot')
@patch('UI_helper.selectManyDataFilesToPlot')
@patch('UI_helper.selectDataFileToPlot')
@patch('UI_helper.selectDataToPlot')

def test(
    mock_selectDataToPlot,
    mock_selectDataFileToPlot,
    mock_selectManyDataFilesToPlot,
    mock_selectObjectToPlot,
    mock_selectFunctionsToRun,
    ):
    mock_selectDataToPlot.return_value = data_name

    mock_selectDataFileToPlot.return_value = 0
    mock_selectManyDataFilesToPlot.return_value = [0,1, int(n_max/2), int(n_max/2)+1, n_max-1, n_max]
    mock_selectObjectToPlot.return_value = (2, 'b')
    # mock_selectFunctionsToRun.return_value = 'all'

    failed_plots = []
    for f in auto_plotters:
        try:
            eval('{}.{}()'.format(f.__module__, f.__name__))
        except:
            failed_plots.append('{}.{}'.format(f.__module__, f.__name__))
    print('failed plotters:')
    for f in failed_plots:
        print(f)

if __name__ == '__main__':
    test()