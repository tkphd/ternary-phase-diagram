# -*- coding:utf-8 -*-

from ctypes import CDLL, c_double
paraboloids = CDLL("./paraboloids.so")

# Control points

xe1A = paraboloids.xe1A
xe2A = paraboloids.xe2A
xe1B = paraboloids.xe1B
xe2B = paraboloids.xe2B
xe1C = paraboloids.xe1C
xe2C = paraboloids.xe2C

xe1A.restype = c_double
xe2A.restype = c_double
xe1B.restype = c_double
xe2B.restype = c_double
xe1C.restype = c_double
xe2C.restype = c_double

# Free energies

GA = paraboloids.GA
GA.argtypes = [c_double, c_double]
GA.restype = c_double

GB = paraboloids.GB
GB.argtypes = [c_double, c_double]
GB.restype = c_double

GC = paraboloids.GC
GC.argtypes = [c_double, c_double]
GC.restype = c_double

# First derivatives

dGAdx1 = paraboloids.dGAdx1
dGAdx1.argtypes = [c_double, c_double]
dGAdx1.restype = c_double

dGAdx2 = paraboloids.dGAdx2
dGAdx2.argtypes = [c_double, c_double]
dGAdx2.restype = c_double

dGBdx1 = paraboloids.dGBdx1
dGBdx1.argtypes = [c_double, c_double]
dGBdx1.restype = c_double

dGBdx2 = paraboloids.dGBdx2
dGBdx2.argtypes = [c_double, c_double]
dGBdx2.restype = c_double

dGCdx1 = paraboloids.dGCdx1
dGCdx1.argtypes = [c_double, c_double]
dGCdx1.restype = c_double

dGCdx2 = paraboloids.dGCdx2
dGCdx2.argtypes = [c_double, c_double]
dGCdx2.restype = c_double

# Second derivatives

d2GAdx11 = paraboloids.d2GAdx11
d2GAdx11.restype = c_double

d2GAdx12 = paraboloids.d2GAdx12
d2GAdx12.restype = c_double

d2GAdx22 = paraboloids.d2GAdx22
d2GAdx22.restype = c_double

d2GBdx11 = paraboloids.d2GBdx11
d2GBdx11.restype = c_double

d2GBdx12 = paraboloids.d2GBdx12
d2GBdx12.restype = c_double

d2GBdx22 = paraboloids.d2GBdx22
d2GBdx22.restype = c_double

d2GCdx11 = paraboloids.d2GCdx11
d2GCdx11.restype = c_double

d2GCdx12 = paraboloids.d2GCdx12
d2GCdx12.restype = c_double

d2GCdx22 = paraboloids.d2GCdx22
d2GCdx22.restype = c_double
