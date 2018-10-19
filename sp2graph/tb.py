"""

:mod:`sp2graph.tb`
===================================

Functions for solving a simple nearest neighbors tight-binding
problem given by the hopping matrix (adjacency matrix).

.. currentmodule:: sp2graph.tb

"""

import numpy as np
from numpy import linalg as LA

__all__ = ['tbBondOrder']


def tbBondOrder(A):
    w, X = LA.eigh(A)

    TB = np.zeros(A.shape, dtype=np.float32)
    for i in range(len(w)):
        if w[i] > 0:
            break
        for j in range(len(A)):
            neig = np.where(A[j, :]==1)[0]
            for k in range(len(neig)):
                idx = neig[k]
                # not sure why the sum is negative...
                TB[j, idx] -= X[j, i]*X[idx, i]
    return TB


def tests():
    """ tests """

if __name__ == '__main__':
    tests()
