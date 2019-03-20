# -*- coding: utf-8 -*-

# Generate phase diagrams
# Usage: python phasediagram.py

from math import ceil, fabs, sqrt
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np
from scipy.optimize import fsolve
from tqdm import tqdm

from pyCinterface import *

density = 251
colors = ['red', 'green', 'blue', 'gray']

# Helper functions to convert compositions into (x,y) coordinates
def simX(x1, x2):
  return x1 + 0.5 * x2

def simY(x2):
  return 0.5 * sqrt(3.) * x2

def euclideanNorm(dx1, dx2):
  return sqrt(dx1**2 + dx2**2)

def boundBy(x, a, b):
  return (a < x) and (x < b)

# Plot phase diagram
pltsize = 10
plt.figure(figsize=(pltsize, 0.5 * sqrt(3.) * pltsize))
plt.title("Ternary Phase Diagram", fontsize=18)
plt.xlabel(r'$x_1$', fontsize=24)
plt.ylabel(r'$x_2$', fontsize=24)
plt.xlim([0, 1])
plt.ylim([0, simY(1)])
# triangle bounding the Gibbs simplex
XS = [0, simX(1,0), simX(0,1), 0]
YS = [0, simY(0),   simY(1),   0]
plt.plot(XS, YS, '-k')

font = FontProperties()
font.set_weight('semibold')
font.set_size(14)
plt.text(simX(0.15, 0.15), simY(0.15), 'A', horizontalalignment='center', verticalalignment='center', fontproperties=font)
plt.text(simX(0.70, 0.15), simY(0.15), 'B', horizontalalignment='center', verticalalignment='center', fontproperties=font)
plt.text(0.5, simY(0.70), 'C', horizontalalignment='center', verticalalignment='center', fontproperties=font)

def ABSolver(x1, x2):
  def system(X):
    x1A, x2A, x1B, x2B = X
    fA = GA(x1A, x2A)
    fB = GB(x1B, x2B)
    dfAdx1 = dGAdx1(x1A, x2A)
    dfAdx2 = dGAdx2(x1A, x2A)
    dfBdx1 = dGBdx1(x1B, x2B)
    dfBdx2 = dGBdx2(x1B, x2B)
    dx1 = x1A - x1B
    dx2 = x2A - x2B
    return [dfAdx1 - dfBdx1,
            dfAdx2 - dfBdx2,
            fA + dfAdx1 * dx1 + dfAdx2 * dx2 - fB,
            (x1 - x1B) * dx2 - dx1 * (x2 - x2B)
           ]

  def jacobian(X):
    x1A, x2A, x1B, x2B = X
    dfAdx1 = dGAdx1(x1A, x2A)
    dfAdx2 = dGAdx2(x1A, x2A)
    dfBdx1 = dGBdx1(x1B, x2B)
    dfBdx2 = dGBdx2(x1B, x2B)
    d2fAdx11 = d2GAdx11()
    d2fAdx12 = d2GAdx12()
    d2fAdx22 = d2GAdx22()
    d2fBdx11 = d2GBdx11()
    d2fBdx12 = d2GBdx12()
    d2fBdx22 = d2GBdx22()
    dx1 = x1A - x1B
    dx2 = x2A - x2B
    return [[ d2fAdx11, d2fAdx12,-d2fBdx11,-d2fBdx12],
            [ d2fAdx12, d2fAdx22,-d2fBdx12,-d2fBdx22],
            [ d2fAdx11 * dx1 + 2*dfAdx1 + d2fAdx12 * dx2,
              d2fAdx12 * dx1 + 2*dfAdx2 + d2fAdx22 * dx2,
              -dfBdx1 - dfAdx1,
              -dfBdx2 - dfAdx2],
            [-x2 + x2B, x1 - x1B, -x2A + x2, -x1 + x1A]
           ]
  # returns the tuple [x1A, x2A, x1B, x2B]
  return fsolve(func=system, x0=[x1, x2, x1, x2], fprime=jacobian)

def ACSolver(x1, x2):
  def system(X):
    x1A, x2A, x1C, x2C = X
    fA = GA(x1A, x2A)
    fC = GC(x1C, x2C)
    dfAdx1 = dGAdx1(x1A, x2A)
    dfAdx2 = dGAdx2(x1A, x2A)
    dfCdx1 = dGCdx1(x1C, x2C)
    dfCdx2 = dGCdx2(x1C, x2C)
    dx1 = x1A - x1C
    dx2 = x2A - x2C
    return [dfAdx1 - dfCdx1,
            dfAdx2 - dfCdx2,
            fA + dfAdx1 * dx1 + dfAdx2 * dx2 - fC,
            (x1 - x1C) * dx2 - dx1 * (x2 - x2C)
           ]

  def jacobian(X):
    x1A, x2A, x1C, x2C = X
    dfAdx1 = dGAdx1(x1A, x2A)
    dfAdx2 = dGAdx2(x1A, x2A)
    dfCdx1 = dGCdx1(x1C, x2C)
    dfCdx2 = dGCdx2(x1C, x2C)
    d2fAdx11 = d2GAdx11()
    d2fAdx12 = d2GAdx12()
    d2fAdx22 = d2GAdx22()
    d2fCdx11 = d2GCdx11()
    d2fCdx12 = d2GCdx12()
    d2fCdx22 = d2GCdx22()
    dx1 = x1A - x1C
    dx2 = x2A - x2C
    return [[ d2fAdx11, d2fAdx12,-d2fCdx11,-d2fCdx12],
            [ d2fAdx12, d2fAdx22,-d2fCdx12,-d2fCdx22],
            [ d2fAdx11 * dx1 + 2*dfAdx1 + d2fAdx12 * dx2,
              d2fAdx12 * dx1 + 2*dfAdx2 + d2fAdx22 * dx2,
              -dfCdx1 - dfAdx1,
              -dfCdx2 - dfAdx2],
            [-x2 + x2C, x1 - x1C, -x2A + x2, -x1 + x1A]
           ]
  # returns the tuple [x1A, x2A, x1C, x2C]
  return fsolve(func=system, x0=[x1, x2, x1, x2], fprime=jacobian)

def BCSolver(x1, x2):
  def system(X):
    x1B, x2B, x1C, x2C = X
    fB = GB(x1B, x2B)
    fC = GC(x1C, x2C)
    dfBdx1 = dGBdx1(x1B, x2B)
    dfBdx2 = dGBdx2(x1B, x2B)
    dfCdx1 = dGCdx1(x1C, x2C)
    dfCdx2 = dGCdx2(x1C, x2C)
    dx1 = x1B - x1C
    dx2 = x2B - x2C
    return [dfBdx1 - dfCdx1,
            dfBdx2 - dfCdx2,
            fB + dfBdx1 * dx1 + dfBdx2 * dx2 - fC,
            (x1 - x1C) * dx2 - dx1 * (x2 - x2C)
           ]

  def jacobian(X):
    x1B, x2B, x1C, x2C = X
    dfBdx1 = dGBdx1(x1B, x2B)
    dfBdx2 = dGBdx2(x1B, x2B)
    dfCdx1 = dGCdx1(x1C, x2C)
    dfCdx2 = dGCdx2(x1C, x2C)
    d2fBdx11 = d2GBdx11()
    d2fBdx12 = d2GBdx12()
    d2fBdx22 = d2GBdx22()
    d2fCdx11 = d2GCdx11()
    d2fCdx12 = d2GCdx12()
    d2fCdx22 = d2GCdx22()
    dx1 = x1B - x1C
    dx2 = x2B - x2C
    return [[ d2fBdx11, d2fBdx12,-d2fCdx11,-d2fCdx12],
            [ d2fBdx12, d2fBdx22,-d2fCdx12,-d2fCdx22],
            [ d2fBdx11 * dx1 + 2*dfBdx1 + d2fBdx12 * dx2,
              d2fBdx12 * dx1 + 2*dfBdx2 + d2fBdx22 * dx2,
              -dfCdx1 - dfBdx1,
              -dfCdx2 - dfBdx2],
            [-x2 + x2C, x1 - x1C, -x2B + x2, -x1 + x1B]
           ]
  # returns the tuple [x1B, x2B, x1C, x2C]
  return fsolve(func=system, x0=[x1, x2, x1, x2], fprime=jacobian)

def ABCSolver(x1, x2):
  def system(X):
    x1A, x2A, x1B, x2B, x1C, x2C = X
    fA = GA(x1A, x2A)
    fB = GB(x1B, x2B)
    fC = GC(x1C, x2C)
    dfAdx1 = dGAdx1(x1A, x2A)
    dfAdx2 = dGAdx2(x1A, x2A)
    dfBdx1 = dGBdx1(x1B, x2B)
    dfBdx2 = dGBdx2(x1B, x2B)
    dfCdx1 = dGCdx1(x1C, x2C)
    dfCdx2 = dGCdx2(x1C, x2C)
    dx1B = x1A - x1B
    dx1C = x1A - x1C
    dx2B = x2A - x2B
    dx2C = x2A - x2C
    return [dfAdx1 - dfBdx1,
            dfAdx1 - dfCdx1,
            dfAdx2 - dfBdx2,
            dfAdx2 - dfCdx2,
            fA + dfAdx1 * dx1B + dfAdx2 * dx2B - fB,
            fA + dfAdx1 * dx1C + dfAdx2 * dx2C - fC
           ]

  def jacobian(X):
    x1A, x2A, x1B, x2B, x1C, x2C = X
    dfAdx1 = dGAdx1(x1A, x2A)
    dfAdx2 = dGAdx2(x1A, x2A)
    dfBdx1 = dGBdx1(x1B, x2B)
    dfBdx2 = dGBdx2(x1B, x2B)
    dfCdx1 = dGCdx1(x1C, x2C)
    dfCdx2 = dGCdx2(x1C, x2C)
    d2fAdx11 = d2GAdx11()
    d2fAdx12 = d2GAdx12()
    d2fAdx22 = d2GAdx22()
    d2fBdx11 = d2GBdx11()
    d2fBdx12 = d2GBdx12()
    d2fBdx22 = d2GBdx22()
    d2fCdx11 = d2GCdx11()
    d2fCdx12 = d2GCdx12()
    d2fCdx22 = d2GCdx22()
    dx1B = x1A - x1B
    dx1C = x1A - x1C
    dx2B = x2A - x2B
    dx2C = x2A - x2C
    return [[ d2fAdx11, d2fAdx12,-d2fBdx11,-d2fBdx12, 0, 0],
            [ d2fAdx11, d2fAdx12,-d2fCdx11,-d2fCdx12, 0, 0],
            [ d2fAdx12, d2fAdx22, 0, 0,-d2fBdx12,-d2fBdx22],
            [ d2fAdx12, d2fAdx22, 0, 0,-d2fCdx12,-d2fCdx22],
            [ d2fAdx11 * dx1B + 2*dfAdx1 + d2fAdx12 * dx2B,
              d2fAdx12 * dx1B + 2*dfAdx2 + d2fAdx22 * dx2B,
              -dfBdx1 - dfAdx1,
              -dfBdx2 - dfAdx2,
              0, 0],
            [ d2fAdx11 * dx1C + 2*dfAdx1 + d2fAdx12 * dx2C,
              d2fAdx12 * dx1C + 2*dfAdx2 + d2fAdx22 * dx2C,
              0, 0,
              -dfCdx1 - dfAdx1,
              -dfCdx2 - dfAdx2]
           ]

  # returns the tuple [x1A, x2A, x1B, x2B]
  return fsolve(func=system, x0=[x1, x2, x1, x2, x1, x2], fprime=jacobian)

pureA = []
pureB = []
pureC = []

tieAB = []
tieAC = []
tieBC = []

triangle = []

for x1test in tqdm(np.linspace(0, 1, density)):
  for x2test in np.linspace(0, 1 - x1test, max(1, ceil((1 - x1test) * density))):
    x1AB, x2AB, x1BA, x2BA = ABSolver(x1test, x2test)
    x1AC, x2AC, x1CA, x2CA = ACSolver(x1test, x2test)
    x1BC, x2BC, x1CB, x2CB = BCSolver(x1test, x2test)

    a = -0.01
    b =  1.01

    ABisPhysical = (boundBy(x1AB, a, b) and boundBy(x2AB, a, b) and
                    boundBy(x1BA, a, b) and boundBy(x2BA, a, b) and
                    boundBy(x1test, min(x1AB, x1BA), max(x1AB, x1BA)) and
                    boundBy(x2test, min(x2AB, x2BA), max(x2AB, x2BA)))
    ACisPhysical = (boundBy(x1AC, a, b) and boundBy(x2AC, a, b) and
                    boundBy(x1CA, a, b) and boundBy(x2CA, a, b) and
                    boundBy(x1test, min(x1AC, x1CA), max(x1AC, x1CA)) and
                    boundBy(x2test, min(x2AC, x2CA), max(x2AC, x2CA)))
    BCisPhysical = (boundBy(x1BC, a, b) and boundBy(x2BC, a, b) and
                    boundBy(x1CB, a, b) and boundBy(x2CB, a, b) and
                    boundBy(x1test, min(x1BC, x1CB), max(x1BC, x1CB)) and
                    boundBy(x2test, min(x2BC, x2CB), max(x2BC, x2CB)))

    # There can be only one three-phase coexistence region.

    if len(triangle) < 1:
      x1ABC, x2ABC, x1BAC, x2BAC, x1CAB, x2CAB = ABCSolver(x1test, x2test)

      ABCisPhysical = (boundBy(x1ABC, a, b) and boundBy(x2ABC, a, b) and
                       boundBy(x1BAC, a, b) and boundBy(x2BAC, a, b) and
                       boundBy(x1CAB, a, b) and boundBy(x2CAB, a, b) and
                       boundBy(x1test, min((x1ABC, x1BAC, x1CAB)), max((x1ABC, x1BAC, x1CAB))) and
                       boundBy(x2test, min((x2ABC, x2BAC, x2CAB)), max((x2ABC, x2BAC, x2CAB))) and
                       boundBy(x1test, min((xe1A(), xe1B(), xe1C())), max((xe1A(), xe1B(), xe1C()))) and
                       boundBy(x2test, min((xe2A(), xe2B(), xe2C())), max((xe2A(), xe2B(), xe2C()))))

      if ABCisPhysical:
        triX = (simX(x1ABC, x2ABC), simX(x1BAC, x2BAC), simX(x1CAB, x2CAB), simX(x1ABC, x2ABC))
        triY = (simY(x2ABC),        simY(x2BAC),        simY(x2CAB),        simY(x2ABC))
        triangle.append((triX, triY))

    # Compute system energies

    fA = GA(x1test, x2test)
    fB = GB(x1test, x2test)
    fC = GC(x1test, x2test)

    fAB = 1e6 * d2GAdx11()
    fAC = 1e6 * d2GAdx11()
    fBC = 1e6 * d2GBdx11()

    if ABisPhysical:
      lAB = euclideanNorm(x1BA - x1AB, x2BA - x2AB)
      wA = euclideanNorm(x1BA - x1test, x2BA - x2test) / lAB
      wB = euclideanNorm(x1test - x1AB, x2test - x2AB) / lAB
      fAB = wA * GA(x1AB, x2AB) + wB * GB(x1BA, x2BA)

    if ACisPhysical:
      lAC = euclideanNorm(x1CA - x1AC, x2CA - x2AC)
      wA = euclideanNorm(x1CA - x1test, x2CA - x2test) / lAC
      wC = euclideanNorm(x1test - x1AC, x2test - x2AC) / lAC
      fAC = wA * GA(x1AC, x2AC) + wC * GC(x1CA, x2CA)

    if BCisPhysical:
      lBC = euclideanNorm(x1CB - x1BC, x2CB - x2BC)
      wB = euclideanNorm(x1CB - x1test, x2CB - x2test) / lBC
      wC = euclideanNorm(x1test - x1BC, x2test - x2BC) / lBC
      fBC = wB * GB(x1BC, x2BC) + wC * GC(x1CB, x2CB)

    energies = np.asarray((fAB, fAC, fBC, fA, fB, fC))
    minIdx = np.argmin(energies)

    if minIdx == 0:
      points = (simX(x1AB, x2AB), simY(x2AB),
                simX(x1BA, x2BA), simY(x2BA))
      tieAB.append(points)
    elif minIdx == 1:
      points = (simX(x1AC, x2AC), simY(x2AC),
                simX(x1CA, x2CA), simY(x2CA))
      tieAC.append(points)
    elif minIdx == 2:
      points = (simX(x1BC, x2BC), simY(x2BC),
                simX(x1CB, x2CB), simY(x2CB))
      tieBC.append(points)
    elif minIdx == 3:
      pureA.append((simX(x1test, x2test), simY(x2test)))
    elif minIdx == 4:
      pureB.append((simX(x1test, x2test), simY(x2test)))
    elif minIdx == 5:
      pureC.append((simX(x1test, x2test), simY(x2test)))

for x, y in triangle:
  plt.plot(x, y, color='black')

for x, y in pureA:
  plt.scatter(x, y, c=colors[0], edgecolor=colors[0], s=1)

for x, y in pureB:
  plt.scatter(x, y, c=colors[1], edgecolor=colors[1], s=1)

for x, y in pureC:
  plt.scatter(x, y, c=colors[2], edgecolor=colors[2], s=1)

for xa, ya, xb, yb in tieAB:
  if boundBy(ya, 0, triangle[0][1][0]):
    plt.scatter(xa, ya, c=colors[0], edgecolor=colors[0], s=1)
    plt.scatter(xb, yb, c=colors[1], edgecolor=colors[1], s=1)
    plt.plot([xa, xb], [ya, yb], color="gray", linewidth=0.1)
  else:
    plt.scatter(xa, ya, c=colors[3], edgecolor=colors[3], s=1)
    plt.scatter(xb, yb, c=colors[3], edgecolor=colors[3], s=1)

for xa, ya, xc, yc in tieAC:
  if boundBy(xa, 0, triangle[0][0][0]):
    plt.scatter(xa, ya, c=colors[0], edgecolor=colors[0], s=1)
    plt.scatter(xc, yc, c=colors[2], edgecolor=colors[2], s=1)
    plt.plot([xa, xc], [ya, yc], color="gray", linewidth=0.1)
  else:
    plt.scatter(xa, ya, c=colors[3], edgecolor=colors[3], s=1)
    plt.scatter(xc, yc, c=colors[3], edgecolor=colors[3], s=1)

for xb, yb, xc, yc in tieBC:
  if boundBy(xb, triangle[0][0][1], 1):
    plt.scatter(xb, yb, c=colors[1], edgecolor=colors[1], s=1)
    plt.scatter(xc, yc, c=colors[2], edgecolor=colors[2], s=1)
    plt.plot([xb, xc], [yb, yc], color="gray", linewidth=0.1)
  else:
    plt.scatter(xb, yb, c=colors[3], edgecolor=colors[3], s=1)
    plt.scatter(xc, yc, c=colors[3], edgecolor=colors[3], s=1)

plt.savefig("ternary-diagram.png", dpi=300, bbox_inches="tight")
