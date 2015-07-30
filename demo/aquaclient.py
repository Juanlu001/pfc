# coding: utf-8
"""Client program.

How to collect data:

1. Positions received as a list of points
2. Interpolate function (probably unstructured grid, see scipy.interpolate.griddata)
3. Evaluate function on points coming from the FEniCS mesh
4. Restore the values onto the array of a FEniCS function in the proper order

Hint: http://fenicsproject.org/qa/3975/interpolating-vector-function-from-python-code-to-fenics#a3976

fe.vector()[V.dofmap().dofs()] = f(x, y)

"""
import sys
import zmq

import pprint

# Socket to talk to server
context = zmq.Context()
socket = context.socket(zmq.REQ)
socket.connect("tcp://localhost:5556")

print("Collecting data...")

dt = 0.1  # s
t_aim = 0.0  # s

# Collect all
while True:
    # We first initialize the server
    socket.send_pyobj({'time': t_aim})

    t_aim += dt

    data = socket.recv_pyobj()
    print(data['time'], data.keys())
