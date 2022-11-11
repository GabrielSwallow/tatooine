import os

all_data_dir = '../Data/'

def selectDataToPlot(datasets = None) -> str:
    if not datasets:
        datasets = os.listdir(all_data_dir)
    for d_index, filename in enumerate(datasets):
        print(d_index, ' : ', filename)
    
    dataset_index = int(input("please choose a dataset to plot for \n"))
    return datasets[dataset_index]