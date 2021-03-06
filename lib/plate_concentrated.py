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
* Compute global and local error with respect to mixed biharmonic solution
* Make plots of the above
* Measure performance
* Extend for other load conditions

"""
from __future__ import division

import numpy as np
from numpy import pi, sin
from scipy import interpolate

from numba import njit

from dolfin import cpp
from dolfin.functions import functionspace, expression
from dolfin.fem import interpolation
from dolfin.common import plotting


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
    # Sum x and y series terms
    for mm in range(1, max_m):
        for nn in range(1, max_n):
            # Sum for every point
            for ii in range(max_i):
                for jj in range(max_j):
                    a_mn = a_mn_point(P, a, b, xi, eta, mm, nn)
                    ww[ii, jj] += (a_mn / (mm**2 / a**2 + nn**2 / b**2)**2
                                   * sin(mm * pi * xx[ii, jj] / a)
                                   * sin(nn * pi * yy[ii, jj] / b)
                                   / (pi**4 * D))


class ExactRectangularPlate(expression.Expression):
    """Exact solution of the transverse displacement of a rectangular plate.

    This expression represents the exact Navier solution of a simply supported
    plate.

    Notes
    -----
    The Navier solution is computed internally using jit-compiled functions
    for higher performance. The result is then interpolated using cubic
    splines.

    """
    def __init__(self, mesh, h, E, nu, P, xi, eta, max_m=64, max_n=64):
        """Constructor.

        Parameters
        ----------
        mesh : mesh.RectangleMesh
            Mesh corresponding to the plate geometry.
        h : float
            Thickness of the plate.
        E : float
            Young modulus [Pa].
        nu : float
            Poisson's ratio.
        P : float
            Point load [N].
        xi, eta : float
            Coordinates of the load application point [m].
        max_m, max_n : int, optional
            Maximum number of series terms to sum, default to 64.

        """
        # Geometry
        assert isinstance(mesh, cpp.mesh.RectangleMesh)
        coords = mesh.coordinates()
        a, b = coords.max(axis=0) - coords.min(axis=0)
        x, y = coords.T
        num_x = len(np.unique(x))
        num_y = len(np.unique(y))
        assert num_x * num_y == mesh.num_vertices()

        # Recreate domain
        # x and y must be swapped, see https://github.com/scipy/scipy/issues/3164
        xx = x.reshape(num_y, num_x)
        yy = y.reshape(num_y, num_x)

        # Material
        D = h**3 * E / (12 * (1 - nu**2))

        # Compute displacement
        ww = np.zeros_like(xx)
        plate_displacement(xx, yy, ww, a, b, P, xi, eta, D, max_m, max_n)

        # Interpolate
        self._w = interpolate.RectBivariateSpline(yy[:, 0], xx[0], ww, kx=3, ky=3)

    def eval(self, value, x):
        value[0] = self._w(x[1], x[0])


if __name__ == '__main__':
    # Plate geometry
    a = 1.0  # m
    b = 2.0  # m
    h = 50e-3  # m

    # Material properties
    E = 69e9  # Pa
    nu = 0.35

    # Load
    P = -10e3  # N
    xi = a / 2
    eta = b / 2

    # Flexural rigidity
    D = h**3 * E / (12 * (1 - nu**2))

    # Now comes the FEniCS magic
    mesh = cpp.mesh.RectangleMesh(0, 0, a, b, 20, 20)
    V = functionspace.FunctionSpace(mesh, 'Lagrange', 1)

    u_expr = ExactRectangularPlate(mesh, h, E, nu, P, xi, eta)
    u = interpolation.interpolate(u_expr, V)

    # Print maximum displacement
    w_max = np.abs(u.vector().array()).max()
    print "Maximum displacement = %14.12f mm" % (w_max * 1e3)
    print "alpha = %7.5f" % (w_max / (P * a**2 / D))

    plotting.plot(u, interactive=True)
