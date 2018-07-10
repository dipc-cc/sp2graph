#!/usr/bin/env python

from __future__ import print_function
import sp2graph.geometry as sp2gge
import sp2graph.visual as sp2gvi
import sp2graph.graph as sp2ggr

syscoord = 'geometries/phenanthrene_rot.xyz'

# read geometry, move to xy plane and assign the vertex array V
V = sp2gge.readgeom(syscoord)
print(V, '\n')

# initial visualization
sp2gvi.viewV(V, 'ini.png', 5, 5)

# create a graph representation and print the adjacency matrix on screen
G = sp2ggr.adjacencyG(V)
sp2gvi.printAdj(G)

# aplly band width reduction and print the adjacency matrix on screen
sp2ggr.reduceBandWidth(G, V)
sp2gvi.printAdj(G)

# visualization after ordering
sp2gvi.viewV(V, 'order.png', 5, 5)
