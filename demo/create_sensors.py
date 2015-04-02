# coding: utf-8
"""Test sensors output.

"""
from __future__ import division

from dolfin import cpp
from dolfin.functions import functionspace, function

a = 0.6  # m
b = 0.15  # m
num_elem = 10

mesh = cpp.mesh.RectangleMesh(0, -b/2, a, b/2, num_elem, num_elem)
mesh.coordinates()
V = functionspace.FunctionSpace(mesh, 'CG', 2)

dofmap = V.dofmap()
dof_coords = dofmap.tabulate_all_coordinates(mesh).reshape((-1, mesh.geometry().dim()))  # N x 3


if __name__ == '__main__':
    with open("Sensors.dat", 'w') as fh:
        for point in dof_coords:
            z, y = point
            x = -1.01
            nx = -1.0
            ny = 0.0
            nz = 0.0
            line = "{} {} {} 0.0, " * 2 + "0.0 0.0 0.0 0.0, " * 2 + "998, 0.0, 0.0, 0\n"
            fh.write(line.format(x, y, z, nx, ny, nz))
