import os, sys
from matplotlib import pyplot as plt
import numpy as np

def makes_lines(nb_lines):
    for i in range(nb_lines):
        print i
        #plt.figure()
        plt.plot(np.arange(10), np.ones(10), linewidth = 20)
        remove_axes()
        plt.savefig('images/fig'+str(i)+'.jpg')

def remove_axes():
    f = plt.gca()
    for tick in f.axes.get_xticklines():
        tick.set_visible(False)
    for tick in f.axes.get_yticklines():
        tick.set_visible(False)
    for xlabel in f.axes.get_xticklabels():
        xlabel.set_visible(False)
        xlabel.set_fontsize(0.0)
    for ylabel in f.axes.get_yticklabels():
        ylabel.set_fontsize(0.0)
        ylabel.set_visible(False)

if __name__ == '__main__':
    makes_lines(10)
    plt.show()