# coding: utf-8
"""Navier solution of a simply supported rectangular plate.

* Accelerated using numba in nopython mode.
* Plotted with FEniCS VTK plotting utility.

Author: Juan Luis Cano Rodríguez <juanlu001@gmail.com>

References
----------
* Timoshenko, S., & Woinowsky-Krieger, S. (1959). "Theory of plates and
  shells" (Vol. 2, p. 120). New York: McGraw-hill.
* Efunda.com, 2014, eFunda: Classical Plate Case Study. [online]. 2014.
  [Accessed 24 December 2014]. Available from:
  http://www.efunda.com/formulae/solid_mechanics/plates/casestudy_list.cfm#SSSS

TODO
----
* Don't hardcode array sizes
* Compute global and local error with respect to mixed biharmonic solution
* Make plots of the above
* Study number of points and interpolation degree
* Check several alpha values from Timoshenko
* Measure performance
* Separate in modules
* Extend for other load conditions
* Don't import all FEniCS

"""
from __future__ import division

import numpy as np
from numpy import pi, sin
from scipy import interpolate

from numba import njit

# If no type hinting is provided, it is inferred in the first call
@njit
def a_mn_point(P, a, b, xi, eta, mm, nn):
    """Navier series coefficient for concentrated load.

    """
    return 4 * P * sin(mm * pi * xi / a) * sin(nn * pi * eta / b) / (a * b)


@njit
def plate_displacement(xx, yy, ww, a, b, P, xi, eta, D, max_m, max_n):
    """Plate transverse displacement.

    """
    # Unroll all the loops for glorious numba efficiency
    max_i, max_j = ww.shape
    for mm in range(1, max_m):
        for nn in range(1, max_n):
            for ii in range(max_i):
                for jj in range(max_j):
                    a_mn = a_mn_point(P, a, b, xi, eta, mm, nn)
                    ww[ii, jj] += (a_mn / (mm**2 / a**2 + nn**2 / b**2)**2
                                   * sin(mm * pi * xx[ii, jj] / a)
                                   * sin(nn * pi * yy[ii, jj] / b)
                                   / (pi**4 * D))


if __name__ == '__main__':
    # --- Initial data

    # Plate geometry
    a = 1.0  # m
    b = 1.0  # m
    h = 50e-3  # m

    # Material properties
    E = 69e9  # Pa
    nu = 0.35

    # Series terms
    max_m = 64
    max_n = 64

    # Computation points
    # NOTE: With an odd number of points the center of the place is included in
    # the grid
    NUM_POINTS = 11

    # Load
    P = -10e3  # N
    xi = a / 2
    eta = a / 2

    # Flexural rigidity
    D = h**3 * E / (12 * (1 - nu**2))

    # ---

    # Set up domain
    x = np.linspace(0, a, num=NUM_POINTS)
    y = np.linspace(0, b, num=NUM_POINTS)
    xx, yy = np.meshgrid(x, y)

    # Compute displacement field
    ww = np.zeros_like(xx)
    plate_displacement(xx, yy, ww, a, b, P, xi, eta, D, max_m, max_n)

    # Don't use RegularGridInterpolator as it is intended for arbitrary dimensions
    # and doesn't support quadratic or cubic interpolation
    w = interpolate.RectBivariateSpline(x, y, ww, kx=3, ky=3)

    # Print maximum displacement
    w_max = abs(ww).max()
    print "Displacement at the center = %14.12f mm" % (w(a / 2, b / 2) * 1e3)
    print "Maximum displacement (abs) = %14.12f mm" % (w_max * 1e3)
    print "alpha = %7.5f" % (w_max / (P * a**2 / D))
    print "alpha * P a^2 / D = %6.4f mm" % (0.01160 * P * a**2 / D * 1e3)

    # Now comes the FEniCS magic
    from fenics import *

    mesh = UnitSquareMesh(100, 100)
    V = FunctionSpace(mesh, 'Lagrange', 1)
    u = Function(V)

    # The numbering of the coordinates in the mesh and the DOFs do not
    # necessarily coincide
    # Source: http://fenicsproject.org/qa/3258/manually-setting-values-to-a-function?show=3259#a3259
    # Workaround: http://fenicsproject.org/qa/2715/coordinates-u_nodal_values-using-numerical-source-function?show=2721#a2721
    dof_coordinates = V.dofmap().tabulate_all_coordinates(mesh)
    n = V.dim()
    d = mesh.geometry().dim()
    dof_coordinates.resize((n, d))
    dof_x = dof_coordinates[:, 0]
    dof_y = dof_coordinates[:, 1]
    # Compute array values
    vals = np.empty((101 * 101), dtype=float)
    for ii in range(len(dof_x)):
        vals[ii] = w(dof_x[ii], dof_y[ii])
    u.vector()[:] = vals

    plot(u, interactive=True)
