"""

:mod:`sp2graph.graph`
===================================

This module contains functions for dealing with the graph
representation of carbon sp2 geometries.

.. currentmodule:: sp2graph.graph

"""

import numpy as np
from numpy import linalg as LA
import sp2graph.linalg_utils as sp2lau
import sys

__all__ = ['adjacencyG', 'adjacencySelfIntG', 'checkPeriodic',
           'periodicDirections', 'revCMK', 'allKekules']


def adjacencyG(V, L=None, radius=1.6):
    """
    Returns the adjacency matrix (which indicates if there
    is an edge, 1, or not, 0, between two vertex) from the
    first neighbours analysis from the vertices in `V`. Case
    the lattice vectors are provided in the matrix `L`,
    analyses the neighbors from periodic cells.
    """
    nV = len(V)
    V = np.array(V)
    G = np.zeros(shape=[nV, nV], dtype=np.uint8)

    # nearest neighbours (inside cell)
    for i in range(nV):
        idx = sp2lau.closeV(i, V, radius)
        G[i, idx] = 1

    # check periodicity
    if np.any(L):
        # neighbors from `V+L[0]` and `V-L[0]`
        for i in range(nV):
            idx = sp2lau.closeV(i, V, radius, L[0])
            G[i, idx] = 1
            G[idx, i] = 1
        # neighbors from `V+L[1]` and `V-L[1]`
        for i in range(nV):
            idx = sp2lau.closeV(i, V, radius, L[1])
            G[i, idx] = 1
            G[idx, i] = 1
        # neighbors from `V+(L[0]+L[1])` `V-(L[0]+L[1])`
        for i in range(nV):
            idx = sp2lau.closeV(i, V, radius, L[0]+L[1])
            G[i, idx] = 1
            G[idx, i] = 1
        # neighbors from `V+(L[0]-L[1])` `V-(L[0]-L[1])`
        for i in range(nV):
            idx = sp2lau.closeV(i, V, radius, L[0]-L[1])
            G[i, idx] = 1
            G[idx, i] = 1
    return G


def periodicDirections(V, L, radius=1.6):
    """
    Check the periodicity and return an `n` by 2 array like `pdir`
    that indicates the `n` periodic directions. Note that we assume
    inversion symmetry, so if the system is periodic in `L[0]` then
    we do not need to include/check `-L[0]`.

    Returns
    -------
    If the system is periodic, for each periodic direction stack
       [1, 0] for periodicity in x,
       [0, 1] for periodicity in y,
       [1, 1] for periodicity in both x and y,
       [1, -1] for periodicity in x and -y
    (note that inversion symmetry is assumed).
    If the system is NOT periodic returns None.
    """
    pdir = []
    # check periodicity
    if np.any(L):
        nV = len(V)
        # check neighbors at `V+L[0]`
        for i in range(nV):
            idx = sp2lau.closeV(i, V, radius, L[0])
            if len(idx):
                pdir.append([1, 0])
                break
        # check neighbors at `V+L[1]`
        for i in range(nV):
            idx = sp2lau.closeV(i, V, radius, L[1])
            if len(idx):
                pdir.append([0, 1])
                break
        # check neighbors at `V+L[0]+L[1]`
        for i in range(nV):
            idx = sp2lau.closeV(i, V, radius, L[0]+L[1])
            if len(idx):
                pdir.append([1, 1])
                break
        # check neighbors at `V+L[0]-L[1]`
        for i in range(nV):
            idx = sp2lau.closeV(i, V, radius, L[0]-L[1])
            if len(idx):
                pdir.append([1, -1])
                break
    return pdir


def checkPeriodic(v1, v2, L, radius=1.6):
    """
    Recieve two vertices `v1` and `v2` and check whether
    they are connected via periodic cells. If so, returns
    `v2` displaced by the corresponding lattice vectors.
    To be used when plotting bonds across periodic cells.
    """

    if np.any(L):
        vp = None
        if LA.norm((v2-v1)) < radius:
            if np.any(vp) == None:
                vp = np.vstack((v1, v2))
            else:
                vp = np.vstack((vp, np.vstack((v1, v2))))
        if LA.norm(v2+L[0]-v1) < radius:
            if np.any(vp) == None:
                vp = np.vstack((np.vstack((v1, v2+L[0])),
                                np.vstack((v1-L[0], v2))))
            else:
                vp = np.vstack((vp,
                                np.vstack((np.vstack((v1, v2+L[0])),
                                           np.vstack((v1-L[0], v2))))))
        if LA.norm(v2-L[0]-v1) < radius:
            if np.any(vp) == None:
                vp = np.vstack((np.vstack((v1, v2-L[0])),
                                np.vstack((v1+L[0], v2))))
            else:
                vp = np.vstack((vp,
                                np.vstack((np.vstack((v1, v2-L[0])),
                                           np.vstack((v1+L[0], v2))))))
        if LA.norm(v2+L[1]-v1) < radius:
            if np.any(vp) == None:
                vp = np.vstack((np.vstack((v1, v2+L[1])),
                                np.vstack((v1-L[1], v2))))
            else:
                vp = np.vstack((vp,
                                np.vstack((np.vstack((v1, v2+L[1])),
                                           np.vstack((v1-L[1], v2))))))
        if LA.norm(v2-L[1]-v1) < radius:
            if np.any(vp) == None:
                vp = np.vstack((np.vstack((v1, v2-L[1])),
                                np.vstack((v1+L[1], v2))))
            else:
                vp = np.vstack((vp,
                                np.vstack((np.vstack((v1, v2-L[1])),
                                           np.vstack((v1+L[1], v2))))))
        if LA.norm(v2+L[0]+L[1]-v1) < radius:
            if np.any(vp) == None:
                vp = np.vstack((np.vstack((v1, v2+L[0]+L[1])),
                                np.vstack((v1-L[0]-L[1], v2))))
            else:
                vp = np.vstack((vp,
                                np.vstack((np.vstack((v1, v2+L[0]+L[1])),
                                           np.vstack((v1-L[0]-L[1], v2))))))
        if LA.norm(v2-L[0]-L[1]-v1) < radius:
            if np.any(vp) == None:
                vp = np.vstack((np.vstack((v1, v2-L[0]-L[1])),
                                np.vstack((v1+L[0]+L[1], v2))))
            else:
                vp = np.vstack((vp,
                                np.vstack((np.vstack((v1, v2-L[0]-L[1])),
                                           np.vstack((v1+L[0]+L[1], v2))))))
        if LA.norm(v2+L[0]-L[1]-v1) < radius:
            if np.any(vp) == None:
                vp = np.vstack((np.vstack((v1, v2+L[0]-L[1])),
                                np.vstack((v1-L[0]+L[1], v2))))
            else:
                vp = np.vstack((vp,
                                np.vstack((np.vstack((v1, v2+L[0]-L[1])),
                                           np.vstack((v1-L[0]+L[1], v2))))))
        if LA.norm(v2-L[0]+L[1]-v1) < radius:
            if np.any(vp) == None:
                vp = np.vstack((np.vstack((v1, v2-L[0]+L[1])),
                                np.vstack((v1+L[0]-L[1], v2))))
            else:
                vp = np.vstack((vp,
                                np.vstack((np.vstack((v1, v2-L[0]+L[1])),
                                           np.vstack((v1+L[0]-L[1], v2))))))
        return vp
    else:
        return np.vstack((v1, v2))


def adjacencySelfIntG(V):
    """
    Returns the adjacency matrix (which indicates if there
    is an edge, 1, or not, 0, between two vertex) from the
    first neighbours analysis from the vertices in V and
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


def reorderResult(R):
    """
    Order a double-bonds list to allow for straight
    comparison among other(s) double-bonds lists.
    """
    DB = np.empty(shape=[0, 0], dtype=np.uint8)
    for i in range(0, len(R)-1, 2):
        if R[i] > R[i+1]:
            foo = R[i]
            R[i] = R[i+1]
            R[i+1] = foo
        if np.all(DB==0):
            DB = [[R[i], R[i+1]]]
        else:
            DB = np.append(DB, [[R[i], R[i+1]]], axis=0)
    DB = DB[DB[:, 0].argsort()]
    return DB


def insertResult(idb, gdb):
    """
    Insert a double-bonds list to the final
    list if not already there.
    """
    for i in range(len(gdb)):
        if np.array_equal(gdb[i], idb):
            return gdb
    gdb = np.append(gdb, [idb], axis=0)
    return gdb


def bfKekules(G, R, Q, DB, rad=None):
    """
    Brute force algorithm that returns in DB all possible Kekule
    structures (i.e., edges with double bonds) from a given adjacency
    matrix G. Vertices with radicals can be provided at 'rad'.
    """

    # TO DO: get rid off this global variable
    global gdb

    # Recursion base: when the queue is empty a search has finished
    idx = len(Q)-1
    if idx < 0:
        # a validy solution has the dimension of the graph
        if len(R) == len(G)-rad.size:
            # append solution to DB if not already there
            idb = reorderResult(R)
            if 'gdb' in globals():
                gdb = insertResult(idb, gdb)
            else:
                gdb = [idb]
            DB = gdb
            return DB
        else:
            allV = np.arange(len(G))
            aux = np.where(np.isin(allV, R, invert=True))[0]
            idx = np.where(np.isin(aux, rad, invert=True))[0]
            Q = np.append(Q, aux[idx[0]])
            return bfKekules(G, R, Q, DB, rad)
    elif Q[0] == -1:
        if 'gdb' in globals():
            DB = gdb
        return DB

    # unqueue the last element (LIFO)
    qval = Q[idx]
    Q = np.delete(Q, idx, 0)

    # when the result array 'R' contains an even number of
    # elements, check whether the unqueued vertex in 'qval'
    # is a valid neigbohr from the the last element on 'R'
    if len(R)%2 != 0 and np.isin(qval, np.where(G[R[-1], :]==1)[0]) == False:
        Q = [-1] # this will stop the search on this branch
        return bfKekules(G, R, Q, DB, rad)

    # append unqueued vertex and analyse its neighbors
    R = np.append(R, qval)
    neig = np.where(G[qval, :]==1)[0]
    # remove from neig those already in R
    dup = np.where(np.isin(neig, R))
    neig = np.delete(neig, dup, 0)
    if rad.size:
        # remove from neig vertices with radicals
        dup = np.where(np.isin(neig, rad))
        neig = np.delete(neig, dup, 0)
    n_neig = len(neig)
    if n_neig > 0:
        # remove from Q those that are also in neig
        dup = np.where(np.isin(Q, neig))[0]
        if len(dup) > 0:
            qdup = Q[dup]
            Q = np.delete(Q, dup, 0)
            dup = np.where(np.isin(neig, qdup))

        if n_neig == 1:
            Q = np.append(Q, neig)
            return bfKekules(G, R, Q, DB, rad)
        else: # we have a bifurcation
            cpQ = Q
            for i in range(n_neig):
                Q = np.append(cpQ, np.roll(neig, i+1))
                if i == n_neig-1:
                    return bfKekules(G, R, Q, DB, rad)
                else:
                    foo = bfKekules(G, R, Q, DB, rad)
    return bfKekules(G, R, Q, DB, rad)


def checkDBlist(G, iniV, DB):
    """
    Check consistency of double bond list
    """
    nDB = len(DB)
    DBt = np.transpose(DB)
    if DBt.shape[0] != 2:
        sys.exit('ERROR: constrained double bonds list has to be tuples!')
    if len(DB.shape) == 1:
        idx = DB[0]
        neig = np.where(G[idx, :]==1)[0]
        if np.isin(DB[1], neig) == False:
            sys.exit('ERROR: vertices %d and %d are not neighbors!'\
                     %(DB[0], DB[1]))
    else:
        for i in range(nDB):
            idx = DB[i, 0]
            neig = np.where(G[idx, :]==1)[0]
            if np.isin(DB[i, 1], neig) == False:
                sys.exit('ERROR: vertices %d and %d are not neighbors!'\
                         %(DB[i, 0], DB[i, 1]))
            if np.any(np.isin(DB[i], DB[0:i])):
                sys.exit('ERROR: double bonds at adjacent edges are not allowed!')
    if np.isin(iniV, DB):
        for i in range(len(G)):
            if np.isin(i, DB) == False:
                iniV = i
                print ('WARNING: initial vertex changed to %d.'%(iniV))
                break
    return iniV


def checkRadlist(G, iniV, C, rad):
    """
    Check consistency of radicals list.
    """
    if C.size:
        if np.any(np.isin(rad, C)):
            sys.exit('ERROR: a vertex assigned as radical cannot belong to the constrained double bonds list!')

    # check if atoms assigned as radicals have only 2 neighbors
    if rad.size == 1:
        if np.sum(G[rad]) > 2:
            sys.exit('ERROR: a radical can only be assigned to carbons having at most 2 neighbors!')
    else:
        for i in rad:
            if np.sum(G[i]) > 2:
                sys.exit('ERROR: a radical can only be assigned to carbons having at most 2 neighbors!')
    if np.isin(iniV, rad):
        allV = np.arange(len(G))
        if C.size:
            allC = C.flatten()
            unc = np.where(np.isin(allV, allC, invert=True))[0]
            if np.any(np.isin(unc, rad)):
                nonrad = np.where(np.isin(unc, rad, invert=True))[0]
                iniV = nonrad[0]
                print ('WARNING: initial vertex changed to %d.'%(iniV))
        else:
            nonrad = np.where(np.isin(allV, rad, invert=True))[0]
            iniV = nonrad[0]
            print ('WARNING: initial vertex changed to %d.'%(iniV))
    return iniV


def allKekules(G, iniV, C=None, rad=None):
    """
    Interface for calling recursive function returning
    all Kekule structures
    """

    if 'gdb' in globals():
        del globals()['gdb']

    R = np.empty(shape=[0], dtype=np.uint8)
    Q = np.empty(shape=[0], dtype=np.uint8)
    DB = np.empty(shape=[0, 0], dtype=np.uint8)

    if C: # constrained search
        C = np.array(C, dtype=np.uint8)
        iniV = checkDBlist(G, iniV, C)
        R = C.flatten()
    else:
        C = np.empty(shape=[0, 0], dtype=np.uint8)

    if rad:
        rad = np.array(rad, dtype=np.uint8)
        iniV = checkRadlist(G, iniV, C, rad)
    else:
        rad = np.empty(shape=[0], dtype=np.uint8)

    Q = np.append(Q, iniV) # queue the starting vertex
    return bfKekules(G, R, Q, DB, rad)


def tests():
    """ tests """

if __name__ == '__main__':
    tests()
