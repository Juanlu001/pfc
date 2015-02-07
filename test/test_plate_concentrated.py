# coding: utf-8
from __future__ import division

from numpy.testing import assert_almost_equal

from fenics import set_log_level, ERROR
set_log_level(ERROR)

from dolfin.cpp import mesh as meshes

from plate_concentrated import SquarePlateDisplacement


def test_different_number_of_elements():
    num_x = 10
    num_y = 20
    mesh = meshes.RectangleMesh(0, 0, 1.0, 1.0, num_x, num_y)

    u_expr = SquarePlateDisplacement(mesh, 1.0, 1.0, 0.1, 1.0, 1.0, 1.0)


def test_timoshenko_factors():
    a = 1.0
    h = 50e-3
    E = 69e9
    nu = 0.35
    P = 10e3

    D = h**3 * E / (12 * (1 - nu**2))

    timoshenko_ratios = (1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 3.0)
    timoshenko_factors = (0.01160, 0.01265, 0.01353, 0.01484, 0.01570,
                          0.01620, 0.01651, 0.01690)
    for b_over_a, factor in zip(timoshenko_ratios, timoshenko_factors):
        b = a * b_over_a
        num_x = 4
        num_y = int(num_x * b_over_a) + int(num_x * b_over_a) % 2
        mesh = meshes.RectangleMesh(0, 0, a, b, num_x, num_y)

        u_expr = SquarePlateDisplacement(mesh, h, E, nu, P, a / 2, b / 2)

        w_max = u_expr(a / 2, b / 2)
        alpha = w_max * D / (P * a**2)

        assert_almost_equal(alpha, factor, decimal=4)
