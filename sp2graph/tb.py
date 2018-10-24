"""

:mod:`sp2graph.tb`
===================================

Functions for solving a simple nearest neighbors tight-binding
problem given by the hopping matrix (adjacency matrix).

.. currentmodule:: sp2graph.tb

"""

import numpy as np
from numpy import linalg as LA
import sp2graph.linalg_utils as sp2lau
import sp2graph.graph as sp2ggr

__all__ = ['tbBondOrder']


def tbBondOrder(V, A, L=None, radius=1.6):
    """ Computes the total bond order (BO) between pairs of neighboring atoms.
    Note that the hopping matrix elements in the electronic Hamiltonian are
    negative numbers as they represent the stabilization energy for electrons
    that are allowed to delocalize, i.e., :math:`H_{ij} < 0`.
    Currently the implementation handles periodic systems along a single
    dimension. It's just a temporary implementation and should be changed
    after the introduction of a proper graph class.
    """
    nV = len(V)

    # nearest neighbours (inside cell)
    A0 = sp2ggr.adjacencyG(V)

    if np.any(L):
        rdir = sp2ggr.periodicDirections(V, L)
        BO = np.zeros(shape=[nV, nV], dtype=np.complex64)
        for i, idir in enumerate(rdir):
            if idir:
                # buld adjacency cell matrices
                Ap1 = np.zeros(shape=[nV, nV], dtype=np.complex64)
                Am1 = np.zeros(shape=[nV, nV], dtype=np.complex64)
                for j in range(nV):
                    idx = sp2lau.closeV(j, V, radius, L[0])
                    Ap1[j, idx] = 1.
                    Am1[idx, j] = 1.

                l = LA.norm(L[i])
                nk = 50 # k-sampling
                dk = 1./nk # step size
                ik = dk/2. # initial k-point
                for k in range(nk):
                    Ak = -1.*np.array(A0, dtype=np.complex64)
                    Ak -= Am1*np.exp(2.0j*np.pi*ik/l)
                    Ak -= Ap1*np.exp(-2.0j*np.pi*ik/l)
                    w, X = LA.eigh(Ak)
                    for i in range(len(w)):
                        if w[i] > 0:
                            break
                        for j in range(nV):
                            # Find neighbor indices
                            neig = np.where(A[j, :]==1)[0]
                            for idx in neig:
                                # Add contribution from pi-orbital:
                                ibo = np.conj(X[j, i])*X[idx, i]
                                BO[j, idx] += ibo + np.conj(ibo)
                    ik += dk
                BO /= nk

        # include the sigma bonds
        BO += 1.*A

        if np.any(rdir):
            return BO.real

    # only gamma point
    w, X = LA.eigh(-1.*A0)
    BO = 1.*A0 # Begin with the sigma bonds
    for i in range(len(w)):
        if w[i] > 0:
            break
        for j in range(nV):
            # Find neighbor indices
            neig = np.where(A0[j, :]==1)[0]
            for idx in neig:
                # Add contribution from pi-orbital (factor 2 for spin):
                BO[j, idx] += 2*X[j, i]*X[idx, i]
    return BO


def tests():
    """ tests """

if __name__ == '__main__':
    tests()
