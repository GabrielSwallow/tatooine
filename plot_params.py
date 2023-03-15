"""
Sets and updates the parameters for plots
"""
import matplotlib.pyplot as plt

def square():
    params = {
        "axes.labelsize":20,
        "font.size":20,
        "legend.fontsize":16,
        "xtick.labelsize":16,
        "ytick.labelsize":16,
        "figure.figsize": [9,8],
    }
    plt.rcParams.update(params)

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