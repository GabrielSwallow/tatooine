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
    if not os.path.isfile('{}/Plots'.format(data_parent_dir)):
        os.mkdir('{}/Plots'.format(data_parent_dir))