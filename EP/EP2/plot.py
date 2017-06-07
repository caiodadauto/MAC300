#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from sympy import symbols
from sympy.plotting import plot

tau = [0.01, 0.05, 0.1, 0.2]
N   = [500, 1000]
x   = np.zeros((4, 20))
y   = np.zeros((4, 20))

for n in N:
    i = 0
    figure(n)
    for t in tau:
        with open(str(n) + "graph-" + str(t)) as f:
            data = f.read()
        data = data.split("\n")
        j = 0
        del data[20]
        for row in data:
            x[i][j] = int(row.split(" ")[0]) + 1
            y[i][j] = float(row.split(" ")[1])
            j += 1
        i += 1
    plt.yscale('log')
    plt.plot(x[0], y[0], 'ro', label='$\\tau = 0.01$')
    plt.plot(x[1], y[1], 'gs', label='$\\tau = 0.05$')
    plt.plot(x[2], y[2], 'b^', label='$\\tau = 0.1$')
    plt.plot(x[3], y[3], 'yv', label='$\\tau = 0.2$')
    grid(color = 'black', linestyle = '--', linewidth = '1', alpha = 0.3)
    legend = plt.legend(loc=3, shadow=True)
    frame = legend.get_frame()
    frame.set_facecolor('0.90')

    for label in legend.get_texts():
        label.set_fontsize('large')

    for label in legend.get_lines():
        label.set_linewidth(1.5)

    plt.xlabel('Iteração')
    plt.ylabel('Residuo')
    plt.subplots_adjust(left = 0.15)
    plt.show()
