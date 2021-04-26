#!/usr/bin/env python

import sp2graph.geometry as sp2gge
import sp2graph.visual as sp2gvi
import sp2graph.graph as sp2ggr
from timeit import default_timer as timer
import glob
import io
from contextlib import redirect_stdout

geoms = sorted(glob.glob('./geometries/graphene-flake_C*'))
for g in geoms:
    # read geometry, move to xy plane and assign the vertex array V
    V, aij = sp2gge.readgeom(g)
    nV = len(V)

    # create a graph representation and print the adjacency matrix on screen
    G = sp2ggr.adjacencyG(V)

    # aplly band width reduction and print the adjacency matrix on screen
    sp2ggr.reduceBandWidth(G, V)

    # calculate all possible Kekule structures (in DB)
    start = timer()
    f = io.StringIO()
    with redirect_stdout(f): # redirect stdout
        for i in range(10): # run 10 times
            Kek = sp2ggr.allKekules(G, 0)
    end = timer()
    time = (end-start)/10. # average over 10 runs
    print(nV, ': ', time, ' s')
