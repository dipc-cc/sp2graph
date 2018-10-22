#!/usr/bin/env python

from __future__ import print_function
import sp2graph.geometry as sp2gge
import sp2graph.visual as sp2gvi
import sp2graph.graph as sp2ggr
import sys

#syscoord = 'geometries/phenanthrene_rot.xyz'
syscoord = sys.argv[1]

# read geometry, move to xy plane and assign the vertex array V
V, aij = sp2gge.readgeom(syscoord)
nV = len(V)
print(nV, '\n')
print(V, '\n')

# create a graph representation and print the adjacency matrix on screen
G = sp2ggr.adjacencyG(V)

# initial visualization
if nV < 200:
    sp2gvi.viewV(V, sizex=5, sizey=5)
    sp2gvi.printAdj(G)

# aplly band width reduction and print the adjacency matrix on screen
sp2ggr.reduceBandWidth(G, V)

# visualization after ordering
if nV < 200:
    sp2gvi.printAdj(G)
    sp2gvi.viewV(V, sizex=5, sizey=5)
