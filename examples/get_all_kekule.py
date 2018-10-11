#!/usr/bin/env python

from __future__ import print_function
import sp2graph.geometry as sp2gge
import sp2graph.visual as sp2gvi
import sp2graph.graph as sp2ggr
import numpy as np

syscoord = 'geometries/phenanthrene.xyz'

# read geometry, move to xy plane and assign the vertex array V
V = sp2gge.readgeom(syscoord)
nV = len(V)
print(nV, '\n')
print(V, '\n')

# create a graph representation and print the adjacency matrix on screen
G = sp2ggr.adjacencyG(V)

# initial visualization
if nV < 100:
    sp2gvi.viewV(V, sizex=5, sizey=5)
    sp2gvi.printAdj(G)

# aplly band width reduction and print the adjacency matrix on screen
sp2ggr.reduceBandWidth(G, V)

# visualization after ordering
if nV < 100:
    sp2gvi.printAdj(G)
    sp2gvi.viewV(V, sizex=5, sizey=5)

R = np.empty(shape=[0], dtype=np.uint8)
Q = np.empty(shape=[0], dtype=np.uint8)
DB = np.empty(shape=[0, 0], dtype=np.uint8)

# start searching from the first vertex with lower degree
ini = np.argmin(sp2ggr.degreeV(V), 0)
Q = np.append(Q, ini)

# calculate all possible Kekule structures (in DB)
DB = sp2ggr.allKekules(G, R, Q, DB)
print('DB ', DB)
