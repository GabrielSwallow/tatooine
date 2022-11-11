import os
import numpy as np

def findMaxFileName(out_dir: str) -> int:
    files = os.listdir(out_dir)
    dataFiles = list(filter(lambda f: f[0:4]=='data', files))
    numbers = [int(d[5:9]) for d in dataFiles]
    maxFileNumber = max(numbers)
    return maxFileNumber

def findNumBodies(out_dir: str) -> int:
    _, pid, _, _ = np.loadtxt( 'nbody.out'.format(out_dir), usecols=(0,1,3,4), unpack=True)
    return len(np.unique(pid))

def createPlotsFolderIfAbsent(data_parent_dir: str) -> None:
    os.path.isfile('/Plots'.format(data_parent_dir))