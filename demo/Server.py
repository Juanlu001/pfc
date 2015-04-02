# coding: utf-8
import numpy as np
import os.path as path
import aquagpusph as aqua

import zmq

n = None
nsensors = None

# Send the data
context = zmq.Context()
socket = context.socket(zmq.REP)
socket.bind("tcp://*:5556")

print "Waiting on recv..."
data = socket.recv_pyobj()
t_aim = data['time']

def main():
    global n, nsensors, t_aim
    t = aqua.get("t")
    if n is None:
        n = 0
        nsensors = 0
        iset = 0
        iset = aqua.get("iset")
        for ii in iset:
            if ii == 0:
                n += 1
            elif ii == 1:
                nsensors += 1
            else:
                raise RuntimeError("Too many isets?")

    # Should we send data to FENICS??
    if t >= t_aim:
        id_inverse = aqua.get("id_inverse")
        p = aqua.get("p")
        pos = aqua.get("r")

        psensors = []
        possensors = []
        for ii in id_inverse[n:n + nsensors]:
            psensors.append(p[ii])  # Ordered pressures
            possensors.append(pos[ii])

        socket.send_pyobj({'time': t, 'p': psensors, 'pos': possensors})

        print "Waiting on recv..."
        data = socket.recv_pyobj()
        t_aim = data['time']

    print "SOCKET SENT"

    return True
