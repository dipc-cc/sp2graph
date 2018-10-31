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

    Parameters
    ----------

    V : array_like
      List of vertices
    A : array_like
      Adjacency matrix
    L : array_like, optional
      Lattice vectors
    radius : float, optional
      Range of neighbor interactions in Angstrom

    Returns
    -------
    numpy.array
      Bond order :math:`\mathrm{BO}_{ij}` between all atom pairs `i` and `j`.
    """
    nV = len(V)

    if np.any(L):
        pdir = sp2ggr.periodicDirections(V, L)
        if np.any(pdir):
            BO = tbBOperiodic(V, L, pdir, radius)

            # include the sigma bonds
            BO += 1.*A

            assert np.allclose(BO.imag, 0*BO)
            return BO.real

    # only gamma point
    w, X = LA.eigh(-1.*A)
    BO = 1.*A # Begin with the sigma bonds
    for i in range(len(w)):
        if w[i] > 0:
            break
        for j in range(nV):
            # Find neighbor indices
            neig = np.where(A[j, :]==1)[0]
            for idx in neig:
                # Add contribution from pi-orbital (factor 2 for spin):
                BO[j, idx] += 2*X[j, i]*X[idx, i]
    return BO


def tbBOperiodic(V, L, pdir, radius):
    """
    Eveluate the Bond Order at a periodic direction.
    """

    nV = len(V)
    BOk = np.zeros(shape=[nV, nV], dtype=np.complex64)

    # nearest neighbours (inside cell)
    A0 = sp2ggr.adjacencyG(V)

    # build adjacency cell matrices to neighboring cells
    npdir = len(pdir)
    Ap1 = np.zeros(shape=[nV, nV, npdir], dtype=np.complex64)
    Am1 = np.zeros(shape=[nV, nV, npdir], dtype=np.complex64)
    for ip, ipdir in enumerate(pdir):
        Ldir = np.dot(ipdir, L)
        for j in range(nV):
            idx = sp2lau.closeV(j, V, radius, Ldir)
            Ap1[j, idx, ip] = 1.
            idx = sp2lau.closeV(j, V, radius, -Ldir)
            Am1[j, idx, ip] = 1.

    # half 1BZ, inv. symmetry for pairs (k, -k)
    kx = (0,)
    for ipdir in pdir:
        if ipdir[0] != 0:
            kx = (np.arange(0.0, 0.5, 0.02)+0.01)
            break
    ky = (0,)
    for ipdir in pdir:
        if ipdir[1] != 0:
            ky = (np.arange(0.0, 0.5, 0.02)+0.01)
            break

    for ikx in kx:
        for iky in ky:
            Ak = -1.*np.array(A0, dtype=np.complex64)
            for ip, ipdir in enumerate(pdir):
                pi2Jnk = 2.0j*np.pi*(ipdir[0]*ikx + ipdir[1]*iky)
                Ak -= Am1[:, :, ip]*np.exp(-pi2Jnk)
                Ak -= Ap1[:, :, ip]*np.exp(pi2Jnk)

            w, X = LA.eigh(Ak)
            for i in range(len(w)):
                if w[i] > 0:
                    break
                for j in range(nV):
                    # Find neighbor indices
                    for idx in np.where(A0[j, :]==1)[0]:
                        # within unit cell
                        ibo = np.conj(X[j, i])*X[idx, i]
                        BOk[j, idx] += ibo + np.conj(ibo)

                    for ip, ipdir in enumerate(pdir):
                        pi2Jnk = 2.0j*np.pi*(ipdir[0]*ikx + ipdir[1]*iky)
                        for idx in np.where(Am1[j, :, ip]==1)[0]:
                            # to neg. side neighbor cell
                            ibo = np.conj(X[j, i])*X[idx, i]*np.exp(-pi2Jnk)
                            BOk[j, idx] += ibo + np.conj(ibo)
                        for idx in np.where(Ap1[j, :, ip]==1)[0]:
                            ibo = np.conj(X[j, i])*X[idx, i]*np.exp(pi2Jnk)
                            BOk[j, idx] += ibo + np.conj(ibo)

    return BOk/(len(kx)*len(ky))


def tests():
    """ tests """

if __name__ == '__main__':
    tests()
