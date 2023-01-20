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

def default():
    plt.rcParams.update(plt.rcParamsDefault)