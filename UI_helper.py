import os
from typing import Tuple
import Navigation_helper
import Data_parser_helper
from Global_variables import *

def selectDataToPlot(datasets = None) -> str:
    if not datasets:
        datasets = os.listdir(all_data_dir)
    for d_index, filename in enumerate(datasets):
        print(d_index, ' : ', filename)
    
    dataset_index = int(input("please choose a dataset to plot for \n"))
    data_name = datasets[dataset_index]
    Navigation_helper.createPlotsFolderIfAbsent(data_name)
    return data_name

def selectDataFileToPlot(out_dir: str) -> int:
    max_file_num = Navigation_helper.findMaxFileName(out_dir)
    print(
        '\nmax file number = {} \nSelect which file to plot (or none to plot {})'
        .format(max_file_num, max_file_num)
        )
    data_file_to_plot = input()
    if data_file_to_plot == '':
        data_file_to_plot = max_file_num
    return int(data_file_to_plot)

def selectManyDataFilesToPlot(out_dir: str):
    max_file_num = Navigation_helper.findMaxFileName(out_dir)
    print(
        '\nmax file number = {} \nSelect which files to plot \neg 1,5,10' # \neg 20-25'
        .format(max_file_num, max_file_num)
        )
    data_files_to_plot = input()
    list_of_files = data_files_to_plot.rsplit(',')
    return [int(f) for f in list_of_files]

def selectFunctionsToRun(functions: list) -> int:
    for index, func in enumerate(functions):
        print('{} : {}.{}'.format(index, func.__module__, func.__name__))
    ui_input = input("please choose a function to run \n")
    if ui_input == 'all':
        return 'all'
    func_index = int(ui_input)
    print('{}.{} selected to run \n'.format(
        functions[func_index].__module__,
        functions[func_index].__name__,
    ))
    return func_index

def selectObjectToPlot(out_dir: str) -> Tuple[int, str]:
    print("\nObjects: \n Smaller Stellar Object : 1 \n Specific Planet : 2, ...")
    obj = int(input("please Select Which Object to Plot \n"))

    num_bodies = Data_parser_helper.findNumBodies(out_dir)
    if obj > num_bodies-1:
        print('\nMax object_id = {}. You selected {}. Try again'.format(num_bodies-1, obj))
        obj, obj_des = selectObjectToPlot(out_dir)

    if obj == 0:
        print('\nCannot select obj = 0 for Nbody characteristics. Try again.')
        obj, obj_des = selectObjectToPlot(out_dir)
    if obj == 1:
        obj_des = 'B'
    elif obj == 2:
        obj_des = 'b'
    elif obj == 3:
        obj_des = 'c'
    elif obj == 4:
        obj_des = 'd'
    else:
        print('\nInvalid obj int entered. Try again.')
        selectObjectToPlot()
        return 

    return obj, obj_des

def selectAnimateRange(out_dir: str):
    max_file_num = Navigation_helper.findMaxFileName(out_dir)
    print(
        '\nmax file number: {} \nSelect range of animation \neg 0,50 \nDefault is 0,{}'.format(max_file_num,max_file_num)
    )
    data_range_raw = input()
    if data_range_raw == '':
        data_range_raw = '0,{}'.format(max_file_num)
    data_range = data_range_raw.rsplit(',')
    return [int(f) for f in data_range]
