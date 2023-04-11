"""
Sets and updates the parameters for plots
"""
import matplotlib.pyplot as plt

def nice_params():
    params = {
        "axes.labelsize":20,
        "font.size":20,
        "legend.fontsize":16,
        "xtick.labelsize":16,
        "ytick.labelsize":16,
        "figure.dpi":500,
    }
    plt.rcParams.update(params)

def square():
    nice_params()
    plt.rcParams.update(
        {"figure.figsize": [7,6]} # [9,8]}
    )

def one_by_two_subplot():
    nice_params()
    plt.rcParams.update(
        {"figure.figsize": [18,8]}
    )

def one_by_N_subplots(N: int):
    nice_params()
    horizontal_size = 7 * N
    plt.rcParams.update(
        {"figure.figsize": [horizontal_size,6]}
    )

def two_by_one_subplot():
    nice_params()
    plt.rcParams.update(
        {"figure.figsize": [9,16]}
    )

def poster():
    params = {
        "axes.labelsize":14,
        "font.size":20,
        "legend.fontsize":14,
        "xtick.labelsize":14,
        "ytick.labelsize":14,
        "figure.figsize":[12,10],
        "figure.dpi":500,
    }
    plt.rcParams.update(params)

def TwoD_sigma_viva():
    params = {
        "axes.labelsize":10,
        "font.size":10,
        "legend.fontsize":12,
        "xtick.labelsize":12,
        "ytick.labelsize":12,
        "figure.figsize":[6,5],
        "figure.dpi":500,
    }
    plt.rcParams.update(params)

def default():
    plt.rcParams.update(plt.rcParamsDefault)