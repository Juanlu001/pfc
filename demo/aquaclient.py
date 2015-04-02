# coding: utf-8
#
# Weather update client
# Connects SUB socket to tcp://localhost:5556
# Collects weather updates and finds avg temp in zipcode
#

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
