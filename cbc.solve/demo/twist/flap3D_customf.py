__author__ = "Harish Narayanan"
__copyright__ = "Copyright (C) 2010 Simula Research Laboratory and %s" % __author__
__license__  = "GNU GPL Version 3 or any later version"

from cbc.twist import *
__VERBOSE__ = False


# Custom expression
class P_F(Expression):
    def __init__(self, problem):
        self._magnitude = 0.5
        self._problem = problem
        self._coordinates = problem.mesh().coordinates()

    def eval_cell(self, value, pos, ufc_cell):
        x = pos[0]
        y = pos[1]
        z = pos[2]
        c_index = ufc_cell.index

        elements = self._problem.elements[c_index]
        e_index = None
        for e in elements:
            p = self._coordinates[e]
            if near(x, p[0]) and near(y, p[1]) and near(z, p[2]):
                e_index = e
                break
        if e_index is None:
            raise ValueError('Failure locating point ({},{},{})'.format(
                x, y, z))

        # FIXME COMPUTE TRANSFORMED/DEFORMED POINT
        # FIXME SEND THE DEFORMED POINT TO SPH

        # Time to know what sensor is in SPH
        try:
            sph_index = self._problem.sph_verts.index(e_index)
            if __VERBOSE__:
                print('element {} -> sensor {}'.format(e_index, sph_index))
        except:
            if __VERBOSE__:
                print('element {} ({}, {}, {}) is not in boundary'.format(
                    e_index, x, y, z))
            """
            value[0] = 0.0
            value[1] = 0.0
            value[2] = 0.0
            return
            """
        # FIXME RETRIEVE THE PRESSURE FROM SPH
        # FIXME COMPUTE THE FORCE

        # print('cell {}, vertex {}'.format(c_index, e_index))

        t = self._problem.solver.t - self._problem.solver.dt

        if(t > 1.0):
            value[0] = self._magnitude
        else:
            value[0] = self._magnitude * t
        value[1] = 0.0
        value[2] = 0.0

    def value_shape(self):
        return (3,)

class Obstruction(Hyperelasticity):

    def mesh(self):
        return Mesh('box.xml')

    def end_time(self):
        return 5.0

    def time_step(self):
        return 0.1

    def is_dynamic(self):
        return True

    def neumann_conditions(self):
        # fluid_force = Expression(("magnitude*min(t, 1.0)", "0.0", "0.0"), magnitude=1.5, t=0)
        fluid_force = P_F(problem=self)
        return [fluid_force]

    def neumann_boundaries(self):
        fluid_interface = "on_boundary && x[2] > -1.0 + DOLFIN_EPS"
        return [fluid_interface]

    def dirichlet_values(self):
        fix = Constant((0.0, 0.0, 0.0))
        return [fix]

    def dirichlet_boundaries(self):
        bottom = "x[2] <= -1.0 + DOLFIN_EPS"
        return [bottom]

    def material_model(self):
        mu    = 60
        lmbda = 90
        material = StVenantKirchhoff([mu, lmbda])
        return material

    def reference_density(self):
        return 1.0

    def time_stepping(self):
        return "CG1"

    def solve(self, elements, sph_verts, sph_verts_coords):
        # Build a dictionary to get the vertexes from a cell
        self.elements = elements
        self.sph_verts = sph_verts
        self.sph_verts_coords = sph_verts
        Hyperelasticity.solve(self)

    def __str__(self):
        return "An obstruction being deformed by an ambient flow"

# Setup problem
problem = Obstruction()

# Get the vertices and facets of the boundary
mesh = problem.mesh()
V = FunctionSpace(mesh, 'CG', 1)
bc = DirichletBC(V, 1, DomainBoundary())
u = Function(V)
bc.apply(u.vector())
d2v = dof_to_vertex_map(V)
vertices_on_boundary = d2v[u.vector() == 1.0]
elements = dict((cell.index(), cell.entities(0)) for cell in cells(mesh))
# Discard the fixed vertices on bottom
coordinates = mesh.coordinates()
verts = []
for vert in vertices_on_boundary:
    coord = coordinates[vert]
    if coord[2] > -1.0 + DOLFIN_EPS:
        verts.append(vert)
del d2v
del bc
del u
del V
del mesh
# Now we can send these vertices to SPH in order to setup the sensors
verts_coord = [coordinates[vert] for vert in verts]
# FIXME SENDING DATA TO SPH

# Solve problem
print problem
problem.solve(elements, verts, verts_coord)
