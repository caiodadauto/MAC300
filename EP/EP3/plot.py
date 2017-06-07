#!/usr/bin/env python
# -*- coding: utf-8 -*-
# python -v 3.5

import math
import numpy as np
import matplotlib.pyplot as plt
from pylab import *

def pol_1(x):
    return (0.448039 + 0.495211 * ((x - 3.878)/2.03) + 0.148228 * (((x - 3.878)/2.03)**2))

def pol_2(x):
    return 41.5347 + 70.8383 * ((x - 16.5)/15.5) + 32.4073 * ((x - 16.5)/15.5)**2

E   = ["exemplo1", "exemplo2"]

i = 0
figure(E[0])
with open(E[0]) as f:
    data = f.read()
data = data.split("\n")
n = int(data[0])
del data[n + 2]
del data[0]
del data[0]
x = np.zeros(n)
y = np.zeros(n)
for row in data:
    x[i] = float(row.split(" ")[0])
    y[i] = float(row.split(" ")[1])
    i += 1
plt.plot(x, pol_1(x), 'b')
plt.plot(x, y, 'ro')
grid(color = 'black', linestyle = '--', linewidth = '1', alpha = 0.3)
plt.xlabel('Tempo (s)')
plt.ylabel('Altura (m)')
plt.subplots_adjust(left = 0.15)


i = 0
figure(E[1])
with open(E[1]) as f:
    data = f.read()
data = data.split("\n")
n = int(data[0])
del data[n + 2]
del data[0]
del data[0]
x = np.zeros(n)
y = np.zeros(n)
for row in data:
    x[i] = float(row.split(" ")[0])
    y[i] = float(row.split(" ")[1])
    i += 1
plt.plot(x, pol_2(x), 'b')
plt.plot(x, y, 'ro')
grid(color = 'black', linestyle = '--', linewidth = '1', alpha = 0.3)
plt.xlabel('Tempo ($\\frac{1}{60}$s)')
plt.ylabel('Posição (cm)')
plt.subplots_adjust(left = 0.15)

plt.show()
