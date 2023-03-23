import os
from typing import Tuple
import Navigation_helper
import Data_parser_helper
from Global_variables import *

# from menu import Menu
# import simple_term_menu as stm

def selectDataToPlot(sub_dir_of_data: str = '', datasets = None) -> str:
    if not datasets:
        datasets = os.listdir(all_data_dir + sub_dir_of_data)
    for d_index, filename in enumerate(datasets):
        print(d_index, ' : ', filename)
    
    dataset_index = int(input("please choose a dataset to plot for \n"))
    data_name = sub_dir_of_data + datasets[dataset_index]
    if 'GROUP' in datasets[dataset_index]:
        return selectDataToPlot(data_name+'/')

    Navigation_helper.createPlotsFolderIfAbsent(data_name)
    return data_name

def selectDataFileToPlot(out_dir: str) -> int:
    # TODO: change this to be data)name not out_dir
    max_file_num = Navigation_helper.findMaxFileNumber(out_dir)
    print(
        '\nmax file number = {} \nSelect which file to plot (or none to plot {})'
        .format(max_file_num, max_file_num)
        )
    data_file_to_plot = input()
    if data_file_to_plot == '':
        data_file_to_plot = max_file_num
    return int(data_file_to_plot)

def selectManyDataFilesToPlot(out_dir: str):
    # TODO: change this to be data)name not out_dir
    max_file_num = Navigation_helper.findMaxFileNumber(out_dir)
    print(
        '\nmax file number = {} \nSelect which files to plot \neg 1,5,10' # \neg 20-25'
        .format(max_file_num, max_file_num)
        )
    data_files_to_plot = input()
    list_of_files = data_files_to_plot.rsplit(',')
    return [int(f) for f in list_of_files]

def select_many_data_ids_to_overlay(choose_out_file: bool = True) -> list[data_id]:
    data_ids = []
    add_more_data_str = 'y'
    while add_more_data_str == 'y':
        data_name = selectDataToPlot()
        directories = Navigation_helper.Directories(data_name)
        if choose_out_file: 
            data_file_to_plot = selectDataFileToPlot(directories.out_dir)
        else: 
            data_file_to_plot = 0
        legend_name = define_legend_name()
        data_ids.append(data_id(data_name, data_file_to_plot, legend_name))
        add_more_data_str = input('continue adding data_ids? \ny for yes \n')
    return data_ids

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

def selectObjectToPlot(out_dir: str) -> astrophysical_object:
    # TODO: change this to be data)name not out_dir
    print("\nObjects: \n Smaller Stellar Object : 1 \n Specific Planet : 2, ...")
    obj_id = int(input("please Select Which Object to Plot \n"))

    num_bodies = Data_parser_helper.findNumBodies(out_dir)
    if obj_id > num_bodies-1:
        print('\nMax object_id = {}. You selected {}. Try again'.format(num_bodies-1, obj_id))
        return selectObjectToPlot(out_dir)

    if obj_id == 0:
        print('\nCannot select obj = 0 for Nbody characteristics. Try again.')
        return selectObjectToPlot(out_dir)

    return objects_in_kep47[obj_id]

def selectObjectsToPlot(out_dir: str) -> Tuple[list, list]:
    print("\nObjects: \n Smaller Stellar Object : 1 \n Specific Planet : 2, ...")
    obj = input("please Select Which Objects to Plot \n")
    obj_list_str = obj.rsplit(',')
    obj_list = []
    for i in obj_list_str:
        obj_list.append(int(i))

    num_bodies = Data_parser_helper.findNumBodies(out_dir)
    obj_des_list = []

    for i in obj_list:
        if i > num_bodies-1:
            print('\nMax object_id = {}. You selected {}. Try again'.format(num_bodies-1, obj))
            obj_list, obj_des_list = selectObjectsToPlot(out_dir)
        if i == 0:
            print('\nCannot select obj = 0 for Nbody characteristics. Try again.')
            obj_list, obj_des_list = selectObjectsToPlot(out_dir)
        obj_des_list.append(objects_in_kep47[i])
    
    return obj_list, obj_des_list

def selectAnimateRange(out_dir: str):
    # TODO: change this to be data)name not out_dir
    max_file_num = Navigation_helper.findMaxFileNumber(out_dir)
    print(
        '\nmax file number: {} \nSelect range of animation \neg 0,50 \nDefault is 0,{}'.format(max_file_num,max_file_num)
    )
    data_range_raw = input()
    if data_range_raw == '':
        data_range_raw = '0,{}'.format(max_file_num)
    data_range = data_range_raw.rsplit(',')
    return [int(f) for f in data_range]

def selectPlottingRange(out_dir: str = None):
    # TODO: change this to be data)name not out_dir
    if out_dir:
        max_file_num = Navigation_helper.findMaxFileNumber(out_dir)
        print('\nmax file number: {} \nSelect range of plot \neg 0,50 \nDefault is 0,{}'.format(max_file_num,max_file_num))
    else:
        print('Select range of plot \neg 0,50 \n max file is not known')
    
    data_range_raw = input()
    if data_range_raw == '':
        data_range_raw = '0,{}'.format(max_file_num)
    data_range = data_range_raw.rsplit(',')
    return [int(f) for f in data_range]

def select_averaging_length() -> int:
    print('Define averaging length (ie number of datapoints to average over). Default is 1, ie none')
    raw_input = input() 
    if raw_input == '':
        return 1
    else:
        return int(raw_input)

def name_the_plot() -> str:
    print('Define the name of the plot:')
    return input()

def define_legend_name() -> str:
    print('Define the name of the data for the legend:')
    return input()


    
# def select_global_vars():
#     print('would you like to update global variables? \n(y/n)')
#     inp = input()
#     if inp == 'y':
#         cart_options = [False, True] # 'cart grid'
#         cart_index = 0
#         logsc_options = [False, True] # log scale
#         logsc_index = 0
#         nbody_options = [False, True] # nbody integrater used
#         nbody_index = 1
#         var_options   = ['rho', 'vx1', 'vx2', 'prs']
#         var_index = 0

#         variables = [
#             '[q] close',
#             '[cart] {}'.format(cart_options[cart_index]),
#             '[logsc] {}'.format(logsc_options[logsc_index]),
#             '[nbody] {}'.format(nbody_options[nbody_index]),
#             '[var] {}'.format(var_options[var_index]),
#             ]
        
#         quitting = False
#         while not quitting:
#             main_menu = stm.TerminalMenu(variables)
#             index = main_menu.show() 
#             optionSelected = variables[index]   
#             if optionSelected == '[q] close':
#                 quitting = True
#             elif optionSelected == '[cart] {}'.format(cart_options[cart_index]):
#                 cart_index = cylce_selection_index(cart_options, cart_index)
#                 variables[index] = '[cart] {}'.format(cart_options[cart_index])
#             elif optionSelected == '[logsc] {}'.format(logsc_options[logsc_index]):
#                 logsc_index = cylce_selection_index(logsc_options, logsc_index)
#                 variables[index] = '[logsc] {}'.format(logsc_options[logsc_index])
#             elif optionSelected == '[nbody] {}'.format(nbody_options[nbody_index]):
#                 nbody_index = cylce_selection_index(nbody_options, nbody_index)
#                 variables[index] = '[nbody] {}'.format(nbody_options[nbody_index])
#             elif optionSelected == '[var] {}'.format(var_options[var_index]):
#                 var_index = cylce_selection_index(var_options, var_index)
#                 variables[index] = '[var] {}'.format(var_options[var_index])

#         with open('Global_variables.py', 'r') as file:
#             string_list = file.readlines()
#         string_list[0] = 'cart = {}\n'.format(cart_options[cart_index])
#         string_list[1] = 'logsc = {}\n'.format(logsc_options[logsc_index])
#         string_list[2] = 'nbody = {}\n'.format(nbody_options[nbody_index])
#         string_list[3] = 'var   = "{}"\n'.format(var_options[var_index])

#         with open('Global_variables.py', 'w') as file:
#             new_file_contents = "".join(string_list)
#             file.write(new_file_contents)

#         return [
#             cart_options[cart_index],
#             logsc_options[logsc_index],
#             nbody_options[nbody_index],
#             var_options[var_index],
#         ]

def cylce_selection_index(options: list, current_index: int) -> int:
    length = len(options)
    if current_index == length - 1:
        return 0
    else: 
        return current_index + 1

if __name__ == '__main__':
    selectDataToPlot()
