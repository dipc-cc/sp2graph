"""

:mod:`sp2graph.graph`
===================================

This module contains functions for dealing with the graph
representation of carbon sp2 geometries.

.. currentmodule:: sp2graph.graph

"""

import numpy as np
import sp2graph.linalg_utils as lau

__all__ = ['adjacencyG', 'adjacencySelfIntG', 'revCMK', 'allKekules']


def adjacencyG(V):
    """
    Returns the adjacency matrix (which indicates if there
    is an edge, 1, or not, 0, between two vertex) from the
    first neighbours analylis from the vertices in V
    """
    nV = len(V)
    G = np.zeros(shape=[nV, nV], dtype=np.uint8)
    for i in range(nV):
        # list nearest neighbours
        idx = lau.closeV(i, V, 1.43)
        # connect nearest neighbours
        G[i, idx] = 1
    return G


def adjacencySelfIntG(V):
    """
    Returns the adjacency matrix (which indicates if there
    is an edge, 1, or not, 0, between two vertex) from the
    first neighbours analylis from the vertices in V and
    include "self-interaction" diagonal term (on elements
    with only two neighbours)
    """
    nV = len(V)
    G = adjacencyG(V)
    for i in range(nV):
        # put 1 in the diagonal when necessary
        G[i, i] = 3-len(idx)
    return G


def degreeV(V):
    """
    Returns an array with the vertices degree
    (2 or 3 in sp2 carbon structures)
    """
    nV = len(V)
    G = adjacencyG(V)
    degV = np.zeros(shape=[nV], dtype=np.uint8)
    for i in range(nV):
        degV[i] = np.sum(G[i])
    return degV


def swapV(G, V, i1, i2):
    """ Swap vertices """
    # swap vertices coordinates
    V[[i1, i2], :] = V[[i2, i1], :]
    # swap columns
    G[:, [i1, i2]] = G[:, [i2, i1]]
    # swap rows
    G[[i1, i2], :] = G[[i2, i1], :]


def naiveReduceBandWidth(G, V):
    """ Naive implementation of band reduction algorithm """
    nV = len(V)
    for i in range(nV-3):
        for j in range(nV-1, i+2, -1):
            if G[i, j]:
                idx = i+2
                while G[i, idx]:
                    idx -= 1
                print(idx, j)
                swapV(G, V, idx, j)


def revCMKinsert(G, degV, labell, queue):
    """ Auxilliary function for revCMK algorithm """
    labell = np.append(labell, queue[0])
    # get neighbours
    ngbr = np.where(G[queue[0], :] == 1)[0]
    # order with increasing degree
    ngbrdeg = np.take(degV, ngbr)
    sortidx = np.unravel_index(np.argsort(ngbrdeg), ngbrdeg.shape)
    ngbr = ngbr[sortidx]
    # remove those already in labell
    mask = np.in1d(ngbr, labell, invert=True)
    ngbr = ngbr[mask]
    # insert neighbours in the queue
    queue = np.append(queue, ngbr)

    return labell, queue


def revCMK(G, V):
    """
    Reverse Cuthill-McKee algorithm implementation.
    Returns an array with the labell order.
    """
    nV = len(V)

    # vertices degrees
    degV = degreeV(V)
    # sorted with increasing degree
    sdegV = np.argsort(degV)

    # result and queue arrays
    labell = np.empty(shape=[0], dtype=np.uint8)
    queue = np.empty(shape=[0], dtype=np.uint8)
    while len(labell) < nV:
        # remove elements already in labell
        mask = np.in1d(sdegV, labell, invert=True)
        sdegV = sdegV[mask]
        queue = np.append(queue, sdegV[0])
        labell, queue = revCMKinsert(G, degV, labell, queue)

        while len(queue) > 0:
            if np.isin(queue[0], labell, invert=True, assume_unique=True):
                labell, queue = revCMKinsert(G, degV, labell, queue)
            # remove the unqueued
            queue = np.delete(queue, 0)

    # reverse
    return labell[::-1]


def reduceBandWidth(G, V):
    """
    Apply the labell order given by a band-width
    reduction algorithm.
    """
    labell = revCMK(G, V)
    nV = len(V)
    idx = np.array(range(nV), dtype=np.uint8)
    for i in range(nV):
        if idx[i] != labell[nV-1-i]:
            j = np.where(idx == labell[nV-1-i])[0]
            swapV(G, V, i, int(j))
            idx[i], idx[j] = idx[j], idx[i]


def allKekules(G, R, Q, DB):
    """
    Brute force algorithm that returns in DB all
    possible Kekule structures (edges with double
    bonds) from a given adjacency matrix G.
    """
    idx = len(Q)-1
    if idx < 0:
        # TO DO HERE:
        # append solution to DB if not already there
        print('R ', R)
        return R
    qval = Q[idx]
    Q = np.delete(Q, idx, 0)
    if qval in R:
        return allKekules(G, R, Q, DB)
    else:
        R = np.append(R, qval)
        neig = np.where(G[qval, :]==1)[0]
        neig = np.sort(neig, axis=None)
        # remove from neig those already in R
        dup = np.where(np.isin(neig, R))
        neig = np.delete(neig, dup, 0)
        n_neig = len(neig)
        if n_neig > 0:
            dup = np.where(np.isin(Q, neig))[0]
            if len(dup) > 0:
                # remove from Q those that are also in neig
                qdup = Q[dup]
                Q = np.delete(Q, dup, 0)
                dup = np.where(np.isin(neig, Q))
                neig = np.delete(neig, dup, 0)
                Q = np.append(Q, neig[::-1])
                Q = np.append(Q, qdup)
            else:
                if n_neig > 1: # we have a bifurcation
                    # OBS.: I think we can consider that we
                    # have at most 2 neighbors here!!!
                    for i in range(n_neig):
                        Q = np.append(Q, np.roll(neig, i+1))
                        if i == n_neig-1:
                            return allKekules(G, R, Q, DB)
                        else:
                            foo = allKekules(G, R, Q, DB)
                else:
                    Q = np.append(Q, neig)
                    return allKekules(G, R, Q, DB)
        return allKekules(G, R, Q, DB)


def tests():
    """ tests """

if __name__ == '__main__':
    tests()
