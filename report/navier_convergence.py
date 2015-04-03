# coding: utf-8
"""Tests exact solution (Navier's) convergence.

"""
from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

from fenics import set_log_level, ERROR
set_log_level(ERROR)

from dolfin.cpp import mesh as meshes
from dolfin.functions import functionspace

from plate_concentrated import ExactRectangularPlate

# Plate geometry
a = 1.0  # m
b = 1.0  # m
h = 50e-3  # m

# Material properties
E = 69e9  # Pa
nu = 0.35

# Load
P = -10e3  # N
xi = a / 2
eta = b / 2

mesh = meshes.RectangleMesh(0, 0, a, b, 10, 10)
V = functionspace.FunctionSpace(mesh, 'Lagrange', 1)

num_iterations = 12
max_w_array = np.zeros(num_iterations)
err_array = np.zeros_like(max_w_array)
for ii in range(num_iterations):
    max_num = 1 << ii

    u_expr = ExactRectangularPlate(mesh, h, E, nu, P, xi, eta,
                                   max_m=max_num, max_n=max_num)

    max_w_array[ii] = np.abs(u_expr(xi, eta))
    err_array[ii] = (max_w_array[ii] - max_w_array[ii - 1]) / (max_w_array[ii] or 1.0)

    print "{:>5d} {:14.12f} {:.3e}".format(max_num, max_w_array[ii], err_array[ii])

print "Plotting results"

max_w_final = max_w_array[-1]
idx_min = np.nonzero(np.abs(max_w_final - max_w_array) < 1e-5)[0][0] + 1

plt.loglog(1 << np.arange(num_iterations)[1:], err_array[1:])
plt.axvline(1 << idx_min, linestyle='--', color='k')
plt.xlabel("Number of series terms")
plt.ylabel("Relative error")
plt.grid()
plt.title("Navier plate convergence")
plt.savefig("navier_convergence.png")
