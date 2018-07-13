#!/usr/bin/env python

from __future__ import print_function
import sp2graph.geometry as sp2gge
import sp2graph.visual as sp2gvi
import sp2graph.graph as sp2ggr
import numpy as np
from timeit import default_timer as timer
import glob
import io
from contextlib import redirect_stdout

geoms = sorted(glob.glob('./geometries/graphene-flake_C*'))
for g in geoms:
    # read geometry, move to xy plane and assign the vertex array V
    V = sp2gge.readgeom(g)
    nV = len(V)

    # create a graph representation and print the adjacency matrix on screen
    G = sp2ggr.adjacencyG(V)

    # aplly band width reduction and print the adjacency matrix on screen
    sp2ggr.reduceBandWidth(G, V)

    R = np.empty(shape=[0], dtype=np.uint8)
    Q = np.empty(shape=[0], dtype=np.uint8)
    DB = np.empty(shape=[0, 0], dtype=np.uint8)

    # start searching from the first vertex with lower degree
    ini = np.argmin(sp2ggr.degreeV(V), 0)
    Q = np.append(Q, ini)

    # calculate all possible Kekule structures (in DB)
    start = timer()
    f = io.StringIO()
    with redirect_stdout(f): # redirect stdout
        for i in range(10): # run 10 times
            R = sp2ggr.allKekules(G, R, Q, DB)
    end = timer()
    time = (end-start)/10. # average over 10 runs
    print(nV, ': ', time, ' s')
