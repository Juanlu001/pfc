# coding: utf-8
"""Tests FEniCS solution convergence.

"""
from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

from dolfin import cpp
from dolfin.functions import function, functionspace, constant
from dolfin.fem import bcs, norms, interpolation, assembling, solving

from ufl import nabla_grad, inner, dx

from plate_concentrated import ExactRectangularPlate

# Plate geometry
a = 1.0  # m
b = 1.0  # m
h = 50e-3  # m

# Material properties
E = 69e9  # Pa
nu = 0.35

# Load
u0 = 0.0  # Displacement at the boundary
M0 = 0.0  # Moment at the boundary
f = 0.0  # Pa, Uniform load
P = -10e3  # N, Centered point load
xi = a / 2
eta = b / 2

D = E * h**3 / (12 * (1 - nu**2))

# Boundary definition
boundary = lambda x, on_boundary: on_boundary

num_iterations = 6
errornorm_array = np.zeros(num_iterations)
for ii in range(2, num_iterations + 2):
    # Initialize geometry
    num_elem = 1 << ii
    mesh = cpp.mesh.RectangleMesh(0, 0, a, b, num_elem, num_elem)

    # FEniCS solution
    CG_u = functionspace.FunctionSpace(mesh, 'CG', 2)
    CG_v = functionspace.FunctionSpace(mesh, 'CG', 2)

    W = CG_u * CG_v  # MixedFunctionSpace

    bc1 = bcs.DirichletBC(W.sub(0), constant.Constant(u0), boundary)
    bc2 = bcs.DirichletBC(W.sub(1), constant.Constant(M0), boundary)
    bc = [bc1, bc2]  # Simply supported

    u, v = function.TrialFunctions(W)
    psi, phi = function.TestFunctions(W)

    lhs = D * (inner(nabla_grad(u), nabla_grad(phi)) +
               inner(nabla_grad(v), nabla_grad(psi)) - inner(v, phi)) * dx
    L = constant.Constant(f) * psi * dx

    point = cpp.mesh.Point(xi, eta)
    P_f = cpp.fem.PointSource(W.sub(0), point, P)

    A, b_v = assembling.assemble_system(lhs, L, bc)
    P_f.apply(b_v)

    w = function.Function(W)
    solving.solve(A, w.vector(), b_v)

    u, v = w.split()  # Displacements, Moments

    # Exact solution
    u_e = ExactRectangularPlate(mesh, h, E, nu, P, xi, eta, max_m=256, max_n=256)
    CG_u_e = functionspace.FunctionSpace(mesh, 'CG', 5)
    u_e_V = interpolation.interpolate(u_e, CG_u_e)

    errornorm_array[ii - 2] = norms.errornorm(u_e_V, u, "l2", mesh=mesh)

    print "{:>3d} {:.3e}".format(num_elem, errornorm_array[ii - 2])

print "Plotting results"

plt.loglog(1 << np.arange(2, num_iterations + 2), errornorm_array)
plt.xlabel("Number of mesh elements")
plt.ylabel("$L_2$ error norm")
plt.grid()
plt.title("FEniCS plate convergence")
plt.savefig("fenics_convergence.png")
