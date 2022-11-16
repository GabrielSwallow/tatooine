import os
import numpy as np

def findMaxFileName(out_dir: str) -> int:
    files = os.listdir(out_dir)
    dataFiles = list(filter(lambda f: f[0:4]=='data', files))
    numbers = [int(d[5:9]) for d in dataFiles]
    maxFileNumber = max(numbers)
    return maxFileNumber

def createPlotsFolderIfAbsent(data_parent_dir: str) -> None:
    os.path.isfile('/Plots'.format(data_parent_dir))
    # TODO: not really necessary, but might be nice