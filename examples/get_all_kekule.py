#!/usr/bin/env python

from __future__ import print_function
import sp2graph.geometry as sp2gge
import sp2graph.visual as sp2gvi
import sp2graph.graph as sp2ggr
import numpy as np

syscoord = 'geometries/naphthalene.xyz'

# read geometry, move to xy plane and assign the vertex array V
V = sp2gge.readgeom(syscoord)
print(V, '\n')

# initial visualization
sp2gvi.viewV(V, sizex=5, sizey=5)

# create a graph representation and print the adjacency matrix on screen
G = sp2ggr.adjacencyG(V)
sp2gvi.printAdj(G)

# aplly band width reduction and print the adjacency matrix on screen
sp2ggr.reduceBandWidth(G, V)
sp2gvi.printAdj(G)

# visualization after ordering
sp2gvi.viewV(V, sizex=5, sizey=5)

R = np.empty(shape=[0], dtype=np.uint8)
Q = np.empty(shape=[0], dtype=np.uint8)
DB = np.empty(shape=[0, 0], dtype=np.uint8)

# start searching from the first vertex with lower degree
ini = np.argmin(sp2ggr.degreeV(V), 0)
Q = np.append(Q, ini)

# calculate all possible Kekule structures (in DB)
R = sp2ggr.allKekules(G, R, Q, DB)
