import os
import Navigation_helper

all_data_dir = '../Data/'

def selectDataToPlot(datasets = None) -> str:
    if not datasets:
        datasets = os.listdir(all_data_dir)
    for d_index, filename in enumerate(datasets):
        print(d_index, ' : ', filename)
    
    dataset_index = int(input("please choose a dataset to plot for \n"))
    return datasets[dataset_index]

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

def selectManyDataFilesToPlot(out_dir: str) -> list[int]:
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
    func_index = int(input("please choose a function to run \n"))
    print('{}.{} selected to run \n'.format(
        functions[func_index].__module__,
        functions[func_index].__name__,
    ))
    return func_index