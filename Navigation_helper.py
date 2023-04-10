import os
import numpy as np
from Global_variables import *
import sys

def findMaxFileNumber(out_dir: str) -> int:
    #TODO: change this to take in data_name
    files = os.listdir(out_dir)
    dataFiles = list(filter(lambda f: f[0:4]=='data', files))
    numbers = [int(d[5:9]) for d in dataFiles]
    maxFileNumber = max(numbers)
    return maxFileNumber

def createPlotsFolderIfAbsent(data_name: str) -> None:
    data_parent_dir = all_data_dir + data_name
    if not os.path.isdir('{}/Plots'.format(data_parent_dir)):
        os.mkdir('{}/Plots'.format(data_parent_dir))

def createDataFolderIfAbsent(data_name: str) -> None:
    data_parent_dir = all_data_dir + data_name
    if not os.path.isdir('{}/data'.format(data_parent_dir)):
        os.mkdir('{}/data'.format(data_parent_dir))

class Directories():
    def __init__(self, data_name, save_plots_local_to_data: bool = False):
        self.all_data_dir = all_data_dir
        self.data_name = data_name
        self.data_parent_dir = all_data_dir + data_name
        if 'plots_for_report.py' in sys.argv[0] and not save_plots_local_to_data:
            self.plots_dir = global_plots_dir
        else: 
            self.plots_dir = self.data_parent_dir + '/Plots/'
        self.out_dir = self.data_parent_dir + '/out'
        self.extra_data_dir = self.data_parent_dir + '/data/'

        self.pluto_log_filename = self.out_dir + '/pluto.log'
        self.nbody_elements_filename = self.out_dir + '/nbody_orbital_elements.out'
        self.nbody_elements_data_filename = self.out_dir + '/nbody_orbital_elements.dat'
        self.nbody_elements_coordinates_filename = self.out_dir + '/nbody_coordinates.dat'
        self.averages = self.out_dir + '/averages.dat'

        self.planet_ini = self.data_parent_dir + '/planet.ini'
    
    def set_data_name(self, new_data_name: str):
        self.data_name = new_data_name