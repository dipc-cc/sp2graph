#!/usr/bin/env python

from __future__ import print_function
import sp2graph.geometry as sp2gge
import sp2graph.visual as sp2gvi
import sp2graph.graph as sp2ggr
import sys

#syscoord = 'geometries/anthracene.xyz'
syscoord = sys.argv[1]

# read geometry, move to xy plane and assign the vertex
# array `V` and lattice vectors in `L` (if provided).
V, L = sp2gge.readgeom(syscoord)
nV = len(V)

# create a graph representation
G = sp2ggr.adjacencyG(V, L)
#sp2ggr.reduceBandWidth(G, V)

# initial visualization of vertices and adjacency matrix
if nV < 100:
    print(nV, '\n')
    print(V, '\n')
    sp2gvi.printAdj(G)

sp2gvi.viewV(V, sizex=5, sizey=5, figname='test_cycles.pdf')
#cycles = sp2ggr.basisCyclesPaton(G)
cycles = sp2ggr.basisCyclesNetworkX(G)

print(cycles)
