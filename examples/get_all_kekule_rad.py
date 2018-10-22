#!/usr/bin/env python

from __future__ import print_function
import sp2graph.geometry as sp2gge
import sp2graph.visual as sp2gvi
import sp2graph.graph as sp2ggr
import sys

#syscoord = 'geometries/triangulene.xyz'
syscoord = sys.argv[1]

# read geometry, move to xy plane and assign the vertex
# array `V` and lattice vectors `aij` (if provided).
V, aij = sp2gge.readgeom(syscoord)
nV = len(V)

# create a graph representation
G = sp2ggr.adjacencyG(V)

# initial visualization of vertices and adjacency matrix
if nV < 10:
    print(nV, '\n')
    print(V, '\n')
    #sp2gvi.viewV(V, sizex=5, sizey=5)
    sp2gvi.printAdj(G)

# apply band width reduction
sp2ggr.reduceBandWidth(G, V)

# vertices and adjacency matrix visualization
if nV < 10:
    #sp2gvi.viewV(V, sizex=5, sizey=5)
    sp2gvi.printAdj(G)

# calculate all possible Kekule structures that contains
# double bonds between the list of edges (tuples) in `C`
C = ((0, 1), (3, 5), (2, 4)) # constrained double bonds
#C = None
# and radical at vertex list in `rad`
#rad = (8, 22) # list of radicals (e.g. for meta-aryne)
rad = (7, 20) # list of radicals (e.g. for triangulene)

Kek = sp2ggr.allKekules(G, 0, C=C, rad=rad)

# visualization of all Kekule structures found
sp2gvi.viewKekuleGrid(V, G, Kek, C=C, rad=rad, sizex=10, sizey=6,
                      figname='kekules_rad.pdf')
#for i in range(len(Kek)):
#    sp2gvi.viewKekule(V, G, Kek[i], C=C, rad=rad, sizex=5, sizey=5)

# visualization of the averaged bond order
sp2gvi.viewBondOrderAverage(V, G, Kek, C=C, rad=rad, sizex=7, sizey=5)
