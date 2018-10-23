#!/usr/bin/env python

from __future__ import print_function
import sp2graph.geometry as sp2gge
import sp2graph.visual as sp2gvi
import sp2graph.graph as sp2ggr
import sys

#syscoord = 'geometries/triangulene.xyz'
syscoord = sys.argv[1]

# read geometry, move to xy plane and assign the vertex
# array `V` and lattice vectors in `L` (if provided).
V, L = sp2gge.readgeom(syscoord)
nV = len(V)

# create a graph representation
G = sp2ggr.adjacencyG(V, L)

# initial visualization of vertices and adjacency matrix
if nV < 10:
    print(nV, '\n')
    print(V, '\n')
    #sp2gvi.viewV(V, sizex=5, sizey=5)
    sp2gvi.printAdj(G)

# calculate all possible Kekule structures that contains
# double bonds between the list of edges (tuples) in `C`
C = ((0, 1)) # constrained double bonds
#C = None
# and radical at vertex list in `rad`
#rad = (25, 4) # list of radicals (e.g. for meta-aryne)
rad = (14, 18) # list of radicals (e.g. for triangulene)

Kek = sp2ggr.allKekules(G, 0, C=C, rad=rad)

# visualization of all Kekule structures found
sp2gvi.viewKekuleGrid(V, G, Kek, C=C, rad=rad, sizex=10, sizey=6,
                      figname='kekules_rad.pdf')
#for i in range(len(Kek)):
#    sp2gvi.viewKekule(V, G, Kek[i], C=C, rad=rad, sizex=5, sizey=5)

# visualization of the Pauling bond order (show constrained
# bonds as usual Kekule representation and radicals as dots)
sp2gvi.viewBondOrderAverage(V, G, Kek, C=C, rad=rad, sizex=7, sizey=5,
                            figname='Pbo_rad.pdf')
