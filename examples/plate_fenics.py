# coding: utf-8
from __future__ import division

from dolfin import cpp
from dolfin.cpp.mesh import Point, RectangleMesh
from dolfin.functions import function, functionspace, constant
from dolfin.fem import bcs, interpolation, assembling, solving
from dolfin.common import plotting

from ufl import nabla_grad, inner, dx

from lib.plate_concentrated import ExactRectangularPlate

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

num_elem = 16

# Boundary definition
boundary = lambda x, on_boundary: on_boundary

# Initialize geometry
mesh = RectangleMesh(Point(0, 0), Point(a, b), num_elem, num_elem)

# FEniCS solution
CG_u = functionspace.FunctionSpace(mesh, 'CG', 2)
CG_v = functionspace.FunctionSpace(mesh, 'CG', 2)

W = CG_u * CG_v  # MixedFunctionSpace

bc1 = bcs.DirichletBC(W.sub(0), constant.Constant(u0), boundary)
bc2 = bcs.DirichletBC(W.sub(1), constant.Constant(M0), boundary)
#bc = [bc1, bc2]  # Simply supported
bc = [bc1]  # Clamped

u, v = function.TrialFunctions(W)
psi, phi = function.TestFunctions(W)

lhs = D * (inner(nabla_grad(u), nabla_grad(phi)) +
           inner(nabla_grad(v), nabla_grad(psi)) - inner(v, phi)) * dx
L = constant.Constant(f) * psi * dx

point = cpp.mesh.Point(xi, eta)
P_f = cpp.fem.PointSource(W.sub(0), point, P)

print("Computing numerical solution...")
A, b_v = assembling.assemble_system(lhs, L, bc)
P_f.apply(b_v)

w = function.Function(W)
solving.solve(A, w.vector(), b_v)

u, v = w.split()  # Displacements, Moments

plotting.plot(u, title="Displacements")
plotting.plot(v, title="Moments")

# Exact solution
u_e = ExactRectangularPlate(mesh, h, E, nu, P, xi, eta, max_m=256, max_n=256)
CG_u_e = functionspace.FunctionSpace(mesh, 'CG', 5)
u_e_V = interpolation.interpolate(u_e, CG_u_e)

# FIXME: UFC is not correctly located
# It is not the first time I encounter this bug, but this is an important one
#plotting.plot(u - u_e_V)

cpp.io.interactive()
