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
nV = len(V)
print(nV, '\n')
print(V, '\n')

# create a graph representation and print the adjacency matrix on screen
G = sp2ggr.adjacencyG(V)

# initial visualization
if nV < 100:
    #sp2gvi.viewV(V, sizex=5, sizey=5)
    sp2gvi.printAdj(G)

# apply band width reduction and print the adjacency matrix on screen
sp2ggr.reduceBandWidth(G, V)

# visualization after ordering
if nV < 100:
    #sp2gvi.viewV(V, sizex=5, sizey=5)
    sp2gvi.printAdj(G)

# calculate all possible Kekule structures (in Kek)
Kek = sp2ggr.allKekules(G, 0)

# visualization of all Kekule structures found
sp2gvi.viewKekuleGrid(V, G, Kek, sizex=10, sizey=6, figname='kekules.pdf')
#for i in range(len(Kek)):
#    sp2gvi.viewKekule(V, G, Kek[i], sizex=5, sizey=5, annotate=True)

# visualization of the averaged bond order
sp2gvi.viewBondOrderAverage(V, G, Kek, sizex=7, sizey=5, figname='bo.pdf', annotate=True)


BO = sp2gtb.tbBondOrder(G)
sp2gvi.viewTBBondOrder(V, BO, sizex=7, sizey=5, figname='bo2.pdf', annotate=True)
