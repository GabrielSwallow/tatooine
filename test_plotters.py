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
import Data_parser_helper

from unittest.mock import patch, create_autospec

disk_properties_plotters = [
    # TwoD_sigma.plot_one, 
    TwoD_sigma.plot_many, 
    # TwoD_sigma.animate,
    # OneD_profile.plot_one, 
    OneD_profile.plot_many, 
    # OneD_profile.animate,
    # Nbody_Characteristics.plot_one_using_dat,
    # Nbody_Characteristics.plot_one_using_out,
    Accretion.plot_one_disk_accretion,
    # Accretion.plot_one_planet_accretion,
    # migration.plot_migration_one,
    # Torque.plot_torque_from_averages,
]

planet_properties_plotters = [
    # Nbody_Characteristics.plot_one_using_dat,
    Nbody_Characteristics.plot_one_using_out,
    Accretion.plot_one_planet_accretion,
    migration.plot_migration_one,
    Torque.plot_torque_from_averages,
]

data_name = UI_helper.selectDataToPlot()
directories = Navigation_helper.Directories(data_name)
n_max = Navigation_helper.findMaxFileNumber(directories.out_dir)

@patch('UI_helper.select_averaging_length')
@patch('UI_helper.selectFunctionsToRun')
# @patch('UI_helper.selectObjectToPlot')
@patch('UI_helper.selectManyDataFilesToPlot')
@patch('UI_helper.selectDataFileToPlot')
@patch('UI_helper.selectDataToPlot')
def plot_disk_properties(
    mock_selectDataToPlot,
    mock_selectDataFileToPlot,
    mock_selectManyDataFilesToPlot,
    # mock_selectObjectToPlot,
    mock_selectFunctionsToRun,
    mock_select_averaging_length
    ):
    mock_selectDataToPlot.return_value = data_name
    mock_select_averaging_length.return_value = 5
    mock_selectDataFileToPlot.return_value = 0
    mock_selectManyDataFilesToPlot.return_value = [0,1, int(n_max/2), int(n_max/2)+1, n_max-1, n_max]
    # mock_selectObjectToPlot.return_value = body_to_plot
    # mock_selectFunctionsToRun.return_value = 'all'

    failed_plots = []
    for f in disk_properties_plotters:
        try:
            eval('{}.{}()'.format(f.__module__, f.__name__))
        except:
            failed_plots.append('{}.{}'.format(f.__module__, f.__name__))
    return failed_plots

@patch('UI_helper.select_averaging_length')
@patch('UI_helper.selectFunctionsToRun')
@patch('UI_helper.selectObjectToPlot')
@patch('UI_helper.selectManyDataFilesToPlot')
@patch('UI_helper.selectDataFileToPlot')
@patch('UI_helper.selectDataToPlot')
def plot_planet_properties(
    object_id,
    mock_selectDataToPlot,
    mock_selectDataFileToPlot,
    mock_selectManyDataFilesToPlot,
    mock_selectObjectToPlot,
    mock_selectFunctionsToRun,
    mock_select_averaging_length,
    ):
    mock_selectDataToPlot.return_value = data_name
    mock_select_averaging_length.return_value = 5
    mock_selectDataFileToPlot.return_value = 0
    mock_selectManyDataFilesToPlot.return_value = [0,1, int(n_max/2), int(n_max/2)+1, n_max-1, n_max]
    mock_selectObjectToPlot.return_value = object_id
    # mock_selectFunctionsToRun.return_value = 'all'

    failed_plots = []
    for f in planet_properties_plotters:
        try:
            eval('{}.{}()'.format(f.__module__, f.__name__))
        except:
            failed_plots.append('{}.{}'.format(f.__module__, f.__name__))
    return failed_plots

if __name__ == '__main__':
    failed_plots = []
    failed_plots.append(plot_disk_properties())
    num_bodies = Data_parser_helper.findNumBodies(directories.out_dir)
    for i in range(2, num_bodies):
        failed_plots.append(plot_planet_properties(objects_ids[i]))
    
    print('failed plotters:')
    for f in failed_plots:
        print(f)