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
    """ Computes the total bond order between pairs of neighboring atoms.
    Note that the hopping matrix elements in the electronic Hamiltonian are
    negative numbers as they represent the stabilization energy for electrons
    that are allowed to delocalize, i.e., :math:`H_{ij} < 0`.
    """
    w, X = LA.eigh(-1.*A)
    TB = 1.*A # Begin with the sigma bonds
    for i in range(len(w)):
        if w[i] > 0:
            break
        for j in range(len(A)):
            # Find neighbor indices
            neig = np.where(A[j, :]==1)[0]
            for idx in neig:
                # Add contribution from pi-orbital (factor 2 for spin):
                TB[j, idx] += 2*X[j, i]*X[idx, i]
    return TB


def tests():
    """ tests """

if __name__ == '__main__':
    tests()
