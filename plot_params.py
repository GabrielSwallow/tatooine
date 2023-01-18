"""
Sets and updates the parameters for plots
"""
import matplotlib.pyplot as plt

def square():
    params = {
        "axes.labelsize":25,
        "font.size":26,
        "legend.fontsize":20,
        "xtick.labelsize":20,
        "ytick.labelsize":20,
        "figure.figsize": [10,10],
    }
    plt.rcParams.update(params)

def default():
    plt.rcParams.update(plt.rcParamsDefault)