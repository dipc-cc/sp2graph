#!/usr/bin/env python

from __future__ import print_function
import sp2graph.geometry as sp2gge
import sp2graph.visual as sp2gvi
import sp2graph.graph as sp2ggr
import sp2graph.tb as sp2gtb
import sys

#syscoord = 'geometries/anthracene.xyz'
syscoord = sys.argv[1]

# read geometry, move to xy plane and assign the vertex array V
V = sp2gge.readgeom(syscoord)

# create a graph representation and print the adjacency matrix on screen
G = sp2ggr.adjacencyG(V)

# apply band width reduction and print the adjacency matrix on screen
sp2ggr.reduceBandWidth(G, V)

# solve the tight-binding problem
TB = sp2gtb.tbBondOrder(G)
sp2gvi.viewTBBondOrder(V, TB, sizex=6, sizey=5)
