"""

:mod:`sp2graph.visual`
===================================

Visualisation of the graph (carbon sp2 structures)

.. currentmodule:: sp2graph.visual

"""

from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
import sp2graph.linalg_utils as lau

__all__ = ['viewV', 'viewKekule', 'printAdj']


def viewV(V, figname=None, sizex=5, sizey=5, dpi=150):
    """ Visualize the vertices (nodes) """
    fig, axs = plt.subplots(figsize=(sizex, sizey))
    nV = len(V)
    for i in range(nV):
        s = 5000/nV
        axs.scatter(V[i, 0], V[i, 1], s=s, c='b', marker=r"$ {} $".format(str(i)), edgecolors='none')
        s = 15000/nV
        axs.scatter(V[i, 0], V[i, 1], s=s, c='b', alpha=.5)
    axs.set_xlim(min(V[:, 0])-2., max(V[:, 0])+2.)
    axs.set_ylim(min(V[:, 1])-2., max(V[:, 1])+2.)
    axs.set_xlabel('x [Ang]')
    axs.set_ylabel('y [Ang]')
    axs.set_aspect('equal')
    if figname:
        plt.savefig(figname, dpi=dpi)
        plt.clf()
        plt.close(fig)
    else:
        plt.show()


def viewKekule(V, A, DB, figname=None, sizex=5, sizey=5, dpi=150):
    """
    Visualize a single Kekule representation with vertices
    coordinates 'V', adjacency matrix 'A' and double-bonds 'DB'.
    """
    fig, axs = plt.subplots(figsize=(sizex, sizey))
    nA = len(A)
    nDB = len(DB)
    for i in range(nA):
        idx = np.transpose(np.nonzero(A[i]))
        for j in range(len(idx)):
            axs.plot((V[i, 0], V[idx[j], 0]),
                     (V[i, 1], V[idx[j], 1]), c='k', ls='-', lw=1.5)
    for i in range(nDB):
        par = lau.parallel(V[DB[i, 0]], V[DB[i, 1]])
        axs.plot((par[0][0], par[1][0]),
                 (par[0][1], par[1][1]), c='r', ls='-', lw=1.5)
    axs.set_xlim(min(V[:, 0])-2., max(V[:, 0])+2.)
    axs.set_ylim(min(V[:, 1])-2., max(V[:, 1])+2.)
    axs.set_xlabel('x [Ang]')
    axs.set_ylabel('y [Ang]')
    axs.set_aspect('equal')
    if figname:
        plt.savefig(figname, dpi=dpi)
        plt.clf()
        plt.close(fig)
    else:
        plt.show()


def printAdj(A):
    """ Print out the adjacency matrix """
    print(' ', end='')
    for i in range(len(A)):
        print('%2d'%i, end='')
    print('\n', A, '\n')


def tests():
    """ tests """
    import numpy as np

    def tviewV():
        V = np.zeros(shape=[1, 2], dtype=np.float32)
        viewV(V)
    tviewV()

if __name__ == '__main__':
    tests()
