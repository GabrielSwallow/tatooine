import os
import numpy as np
from Global_variables import *

def findMaxFileName(out_dir: str) -> int:
    files = os.listdir(out_dir)
    dataFiles = list(filter(lambda f: f[0:4]=='data', files))
    numbers = [int(d[5:9]) for d in dataFiles]
    maxFileNumber = max(numbers)
    return maxFileNumber

def createPlotsFolderIfAbsent(data_name: str) -> None:
    data_parent_dir = all_data_dir + data_name
    if not os.path.isdir('{}/Plots'.format(data_parent_dir)):
        os.mkdir('{}/Plots'.format(data_parent_dir))

class Directories():
    def __init__(self, data_name: str | None) -> None:
        self.all_data_dir = '../Data/'
        self.data_name = data_name
        self.data_parent_dir = all_data_dir + data_name
        self.plots_dir = self.data_parent_dir + '/Plots/'
        self.out_dir = self.data_parent_dir + '/out'

        self.pluto_log_filename = self.out_dir + '/pluto.log'
        self.nbody_elements_filename = self.out_dir + '/nbody_orbital_elements.out'
        self.nbody_elements_data_filename = self.out_dir + '/nbody_orbital_elements.dat'
        self.planet_ini = self.data_parent_dir + '/planet.ini'

    
    def set_data_name(self, new_data_name: str):
        self.data_name = new_data_name