#!/usr/bin/env python

from __future__ import print_function
import sp2graph.geometry as sp2gge
import sp2graph.visual as sp2gvi
import sp2graph.graph as sp2ggr
import sp2graph.linalg_utils as sp2ggla
import numpy as np
import sys

#syscoord = 'geometries/anthracene.xyz'
syscoord = sys.argv[1]

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

# apply band width reduction and print the adjacency matrix on screen
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

# visualization of all Kekule structures found
for i in range(len(DB)):
    sp2gvi.viewKekule(V, G, DB[i], sizex=5, sizey=5)
