#! /usr/bin/env python
#########################################################################
#                                                                       #
#            #    ##   #  #   #                           #             #
#           # #  #  #  #  #  # #                          #             #
#          ##### #  #  #  # #####  ##  ###  #  #  ## ###  ###           #
#          #   # #  #  #  # #   # #  # #  # #  # #   #  # #  #          #
#          #   # #  #  #  # #   # #  # #  # #  #   # #  # #  #          #
#          #   #  ## #  ##  #   #  ### ###   ### ##  ###  #  #          #
#                                    # #             #                  #
#                                  ##  #             #                  #
#                                                                       #
#########################################################################
#
#  This file is part of AQUA-gpusph, a free CFD program based on SPH.
#  Copyright (C) 2012  Jose Luis Cercos Pita <jl.cercos@upm.es>
#
#  AQUA-gpusph is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  AQUA-gpusph is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with AQUA-gpusph.  If not, see <http://www.gnu.org/licenses/>.
#
#########################################################################

import os.path as path
import math

# Input data
# ==========

g = 9.81
hfac = 2.0
cs = 45.0
courant = 0.2
gamma = 7.0
refd = 998.0
alpha = 0.001
delta = 1.0
visc_dyn = 0.000894
# Initial fluid dimensions
h = 0.3
l = 0.6
d = 0.15
# Tank dimensions
H = 0.6
L = 1.61
D = d
# Stimated required number of fluid particles
n = 10000

# Number of sensors (read it from somewhere!!!!)
sensor_data = open("Sensors.dat", 'r').readlines()
n_sensors = len(sensor_data)

# Dimensions and number of particles readjustment
# ===============================================

Vol = l * d * h
dv = Vol / n
dr = dv**(1.0 / 3.0)

nx = int(round(l / dr))
ny = int(round(d / dr))
nz = int(round(h / dr))

hFluid = nz * dr
visc_dyn = max(alpha / 10.0 * refd * hfac * dr * cs, visc_dyn)

Nx = int(round(L / dr))
Ny = int(round(D / dr))
Nz = int(round(H / dr))

prb = cs * cs * refd / gamma

output = open("Fluid.dat", "w")
string = """#############################################################
#                                                           #
#    #    ##   #  #   #                           #         #
#   # #  #  #  #  #  # #                          #         #
#  ##### #  #  #  # #####  ##  ###  #  #  ## ###  ###       #
#  #   # #  #  #  # #   # #  # #  # #  # #   #  # #  #      #
#  #   # #  #  #  # #   # #  # #  # #  #   # #  # #  #      #
#  #   #  ## #  ##  #   #  ### ###   ### ##  ###  #  #      #
#                            # #             #              #
#                          ##  #             #              #
#                                                           #
#############################################################
"""
output.write(string)
print(string)

# Particles generation
# ====================

x0 = 0.0
y0 = - 0.5 * Ny * dr
z0 = 0.0

print('Fluid particles...')
n = 0
Percentage = -1
for i in range(nx):
    if Percentage != (n * 100) / (nx * ny * nz):
        Percentage = (n * 100) / (nx * ny * nz)
        if not Percentage % 10:
            string = '    {}%'.format(Percentage)
            print(string)
    x = x0 + (i + 0.5) * dr
    for j in range(ny):
        y = y0 + (j + 0.5) * dr
        for k in range(nz):
            z = z0 + (k + 0.5) * dr
            imove = 1
            press = refd * g * (hFluid - z)
            dens = pow(press / prb + 1.0, 1.0 / gamma) * refd
            mass = dens * dr**3.0
            string = ("{} {} {} 0.0, " * 4 + "{}, {}, {}, {}\n").format(
                x, y, z,
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0,
                dens,
                0.0,
                mass,
                imove)
            output.write(string)
            n += 1

# Tank generation
# ===============

x0 = -(Nx - nx) * dr
y0 = - 0.5 * Ny * dr
z0 = 0.0
dens = refd
press = 0.0
mass = dr**2.0
print('Tank boundary elements...')
N = 0

print('    Bottom face...')
normal = (0.0, 0.0, -1.0)
z = z0
for i in range(Nx):
    x = x0 + (i + 0.5) * dr
    for j in range(Ny):
        y = y0 + (j + 0.5) * dr
        imove = -3
        string = ("{} {} {} 0.0, " * 4 + "{}, {}, {}, {}\n").format(
            x, y, z,
            normal[0], normal[1], normal[2],
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            dens,
            0.0,
            mass,
            imove)
        output.write(string)
        N += 1

print('    Top face...')
normal = (0.0, 0.0, 1.0)
z = z0 + Nz * dr
for i in range(Nx):
    x = x0 + (i + 0.5) * dr
    for j in range(Ny):
        y = y0 + (j + 0.5) * dr
        imove = -3
        string = ("{} {} {} 0.0, " * 4 + "{}, {}, {}, {}\n").format(
            x, y, z,
            normal[0], normal[1], normal[2],
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            dens,
            0.0,
            mass,
            imove)
        output.write(string)
        N += 1

print('    Front face...')
normal = (-1.0, 0.0, 0.0)
x = x0
for i in range(Ny):
    y = y0 + (i + 0.5) * dr
    for j in range(Nz):
        z = z0 + (j + 0.5) * dr
        imove = -3
        string = ("{} {} {} 0.0, " * 4 + "{}, {}, {}, {}\n").format(
            x, y, z,
            normal[0], normal[1], normal[2],
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            dens,
            0.0,
            mass,
            imove)
        output.write(string)
        N += 1

print('    Back face...')
normal = (1.0, 0.0, 0.0)
x = x0 + Nx * dr
for i in range(Ny):
    y = y0 + (i + 0.5) * dr
    for j in range(Nz):
        z = z0 + (j + 0.5) * dr
        imove = -3
        string = ("{} {} {} 0.0, " * 4 + "{}, {}, {}, {}\n").format(
            x, y, z,
            normal[0], normal[1], normal[2],
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            dens,
            0.0,
            mass,
            imove)
        output.write(string)
        N += 1

print('    Left face...')
normal = (0.0, -1.0, 0.0)
y = y0
for i in range(Nx):
    x = x0 + (i + 0.5) * dr
    for j in range(Nz):
        z = z0 + (j + 0.5) * dr
        imove = -3
        string = ("{} {} {} 0.0, " * 4 + "{}, {}, {}, {}\n").format(
            x, y, z,
            normal[0], normal[1], normal[2],
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            dens,
            0.0,
            mass,
            imove)
        output.write(string)
        N += 1

print('    Right face...')
normal = (0.0, 1.0, 0.0)
y = y0 + Ny * dr
for i in range(Nx):
    x = x0 + (i + 0.5) * dr
    for j in range(Nz):
        z = z0 + (j + 0.5) * dr
        imove = -3
        string = ("{} {} {} 0.0, " * 4 + "{}, {}, {}, {}\n").format(
            x, y, z,
            normal[0], normal[1], normal[2],
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            dens,
            0.0,
            mass,
            imove)
        output.write(string)
        N += 1

output.close()
print('OK')

# XML definition generation
# =========================

templates_path = path.join('./', 'templates')
XML = ('Fluids.xml', 'Main.xml', 'Sensors.xml', 'Settings.xml', 'SPH.xml',
       'Time.xml', 'h_sensor.cl', 'Server.xml')

domain_min = (-(L - l + 0.1), -(0.5 * D + 0.1), -0.1, 0.0)
domain_min = str(domain_min).replace('(', '').replace(')', '')
domain_max = (l + 0.1, 0.5 * D + 0.1, H + 0.1, 0.0)
domain_max = str(domain_max).replace('(', '').replace(')', '')

data = {'DR':str(dr), 'HFAC':str(hfac), 'CS':str(cs), 'COURANT':str(courant),
        'DOMAIN_MIN':domain_min, 'DOMAIN_MAX':domain_max, 'GAMMA':str(gamma),
        'REFD':str(refd), 'VISC_DYN':str(visc_dyn), 'DELTA':str(delta),
        'G':str(g), 'N_SENSORS':str(n_sensors), 'N':str(n + N)}
for fname in XML:
    # Read the template
    f = open(path.join(templates_path, fname), 'r')
    txt = f.read()
    f.close()
    # Replace the data
    for k in data.keys():
        txt = txt.replace('{{' + k + '}}', data[k])
    # Write the file
    f = open(fname, 'w')
    f.write(txt)
    f.close()
