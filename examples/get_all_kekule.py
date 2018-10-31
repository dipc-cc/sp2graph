#!/usr/bin/env python

from __future__ import print_function
import sp2graph.geometry as sp2gge
import sp2graph.visual as sp2gvi
import sp2graph.graph as sp2ggr
import sp2graph.tb as sp2gtb
import sys

#syscoord = 'geometries/anthracene.xyz'
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

# calculate all possible Kekule structures (in Kek)
Kek = sp2ggr.allKekules(G, 0)

# visualization of all Kekule structures found
sp2gvi.viewKekuleGrid(V, G, Kek, sizex=10, sizey=6,
                      figname='kekules.pdf')
#for i in range(len(Kek)):
#    sp2gvi.viewKekule(V, G, Kek[i], sizex=5, sizey=5)

# visualization of the Pauling bond order
sp2gvi.viewBondOrderAverage(V, G, Kek, sizex=7, sizey=5,
                            annotate=True, astyle='white',
                            figname='Pbo.pdf')

# visualization of the Huckel bond order
BO = sp2gtb.tbBondOrder(V)
sp2gvi.viewTBBondOrder(V, BO, sizex=7, sizey=5,
                       annotate=True, astyle='boxed',
                       figname='Hbo.pdf')
