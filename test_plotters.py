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
from analyse import plotters

from unittest.mock import patch, create_autospec


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

    mock_selectDataToPlot.return_value = '1000_Orbits_full'
    mock_selectDataFileToPlot.return_value = 0
    mock_selectManyDataFilesToPlot.return_value = [0,1,2]
    mock_selectObjectToPlot.return_value = 2
    mock_selectFunctionsToRun.return_value = 'all'

    for f in plotters:
        eval('{}.{}()'.format(f.__module__, f.__name__))

if __name__ == '__main__':
    test()