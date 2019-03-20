# -*- coding: utf-8 -*-

# Generate ternary phase diagram
# Usage: python ternary-landscape.py

from math import ceil, sqrt
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from scipy.optimize import fsolve
from tqdm import tqdm

from pyCinterface import *

pltsize = 10
density = 200
ncontour = 70

def simX(x1, x2):
  return x1 + 0.5 * x2

def simY(x2):
  return 0.5 * sqrt(3.) * x2

plt.figure(figsize=(pltsize, 0.5 * sqrt(3.) * pltsize))
plt.axis('off')

x = []
y = []
z = []

for x1test in tqdm(np.linspace(0, 1, density)):
  for x2test in np.linspace(0, 1 - x1test, max(1, ceil((1 - x1test) * density))):
    fA = GA(x1test, x2test)
    fB = GB(x1test, x2test)
    fC = GC(x1test, x2test)

    energies = np.asarray([fA, fB, fC])
    minIdx = np.argmin(energies)

    x.append(simX(x1test, x2test))
    y.append(simY(x2test))
    z.append(energies[minIdx])

fmin = min(z)
fmax = max(z)
print "Raw data spans [{0:2.2e}, {1:2.2e}].".format(fmin, fmax)

x = np.asarray(x)
y = np.asarray(y)
z = np.asarray(z)

levels = np.linspace(0, fmax, ncontour)

plt.tricontourf(x, y, z, levels, cmap=plt.cm.get_cmap('binary'))

plt.savefig("ternary-landscape.png", dpi=300, bbox_inches="tight", transparent=True)
