"""

:mod:`sp2graph.visual`
===================================

Visualisation of the graph (carbon sp2 structures)

.. currentmodule:: sp2graph.visual

"""

from __future__ import print_function
import matplotlib.pyplot as plt

__all__ = ['viewV', 'printAdj']

def viewV(V, figname, sizex=4, sizey=10, dpi=150):
    """ Visualize the vertices (nodes) """
    fig, axs = plt.subplots(figsize=(sizex, sizey))
    for i in range(len(V)):
        axs.scatter(V[i, 0], V[i, 1], s=200, c='b', marker=r"$ {} $".format(str(i)), edgecolors='none')
        axs.scatter(V[i, 0], V[i, 1], s=1000, c='b', alpha=.5)
    axs.set_xlim(min(V[:, 0])-2., max(V[:, 0])+2.)
    axs.set_ylim(min(V[:, 1])-2., max(V[:, 1])+2.)
    axs.set_xlabel('x [Ang]')
    axs.set_ylabel('y [Ang]')
    axs.set_aspect('equal')
    plt.savefig(figname, dpi=dpi);
    plt.clf();
    plt.close(fig)


def printAdj(A):
    """ Print out the adjacency matrix """
    print(' ', end='')
    for i in range(len(A)):
        print('%2d'%i, end='')
    print('')
    print(A)


def tests():
    """ tests """
    import numpy as np

    def tviewV():
        V = np.zeros(shape=[1, 2], dtype=np.float32)
        viewV(V)

if __name__ == '__main__':
    tests()
