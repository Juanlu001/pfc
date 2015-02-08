# coding: utf-8
from dolfin import cpp
from dolfin.functions import function, functionspace, constant
from dolfin.fem import bcs, assembling, solving
from dolfin.common import plotting

from ufl import nabla_grad, inner, dx

# --- Problem data
a = 1.0  # m
b = 2.0  # m
h = 0.050  # m

E = 69e9  # Pa
nu = 0.35

u0 = 0.0  # Displacement at the boundary
M0 = 0.0  # Moment at the boundary
f = 0.0  # Pa, Uniform load
P = -10e3  # N, Centered point load
# ---

num_elem = 20
mesh = cpp.mesh.RectangleMesh(0, 0, a, b, num_elem, num_elem)

CG_u = functionspace.FunctionSpace(mesh, 'CG', 2)
CG_v = functionspace.FunctionSpace(mesh, 'CG', 2)

W = CG_u * CG_v  # MixedFunctionSpace

boundary = lambda x, on_boundary: on_boundary

bc1 = bcs.DirichletBC(W.sub(0), constant.Constant(u0), boundary)
bc2 = bcs.DirichletBC(W.sub(1), constant.Constant(M0), boundary)
bc = [bc1, bc2]  # Simply supported
# bc = [bc1]  # Clamped

u, v = function.TrialFunctions(W)
psi, phi = function.TestFunctions(W)

D = E * h**3 / (12 * (1 - nu**2))

lhs = D * (inner(nabla_grad(u), nabla_grad(phi)) +
           inner(nabla_grad(v), nabla_grad(psi)) - inner(v, phi)) * dx
L = constant.Constant(f) * psi * dx

point = cpp.mesh.Point(a / 2, b / 2)
P_f = cpp.fem.PointSource(W.sub(0), point, P)

A, b = assembling.assemble_system(lhs, L, bc)
P_f.apply(b)

w = function.Function(W)
solving.solve(A, w.vector(), b)

u, v = w.split()

plotting.plot(u, interactive=True)
