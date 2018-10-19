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

__all__ = ['viewV', 'viewKekule', 'viewKekuleGrid',
           'viewBondOrderAverage', 'viewTBBondOrder', 'printAdj']


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


def viewKekule(V, A, DB, figname=None, sizex=5, sizey=5, dpi=150, annotate=False):
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
    if annotate:
        for i in range(nA):
            axs.annotate(i, (V[i, 0], V[i, 1]))
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


def viewKekuleGrid(V, A, DB, figname=None, sizex=5, sizey=5, dpi=150, annotate=False):
    """
    Visualize a single Kekule representation with vertices
    coordinates 'V', adjacency matrix 'A' and double-bonds 'DB'.
    """
    nA = len(A)
    nKek = len(DB)
    nDB = DB.shape[1]

    # define the grid for ploting the figures
    molsize = max(V[:, 0]) - min(V[:, 0]) + 4.
    cols = int(40//molsize)
    rows = nKek/cols + nKek%cols

    fig = plt.figure()
    fig.set_size_inches(sizex, sizey)
    for idb in range(nKek):

        plt.subplot2grid((rows, cols), (idb/cols, idb%cols))
        for i in range(nA):
            idx = np.transpose(np.nonzero(A[i]))
            for j in range(len(idx)):
                plt.plot((V[i, 0], V[idx[j], 0]),
                         (V[i, 1], V[idx[j], 1]), c='k', ls='-', lw=1.5)
        kek = DB[idb]
        for i in range(nDB):
            par = lau.parallel(V[kek[i, 0]], V[kek[i, 1]])
            plt.plot((par[0][0], par[1][0]),
                     (par[0][1], par[1][1]), c='r', ls='-', lw=1.5)
        if annotate:
            for i in range(nA):
                plt.annotate(i, (V[i, 0], V[i, 1]))

        if idb%cols == 0:
            plt.ylabel('y [Ang]')
        else:
            plt.ylabel('')
            plt.tick_params(labelleft='off')
        plt.xlim(min(V[:, 0])-2., max(V[:, 0])+2.)
        plt.ylim(min(V[:, 1])-2., max(V[:, 1])+2.)
        plt.xlabel('x [Ang]')
        plt.gca().set_aspect(1)

    if figname:
        plt.savefig(figname, dpi=dpi)
        plt.clf()
        plt.close(fig)
    else:
        plt.show()


def viewBondOrderAverage(V, A, DB, figname=None, sizex=5, sizey=5, dpi=150, annotate=False):
    """
    Visualize a single Kekule representation with vertices
    coordinates 'V', adjacency matrix 'A' and double-bonds 'DB'.
    """
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    fig, axs = plt.subplots(figsize=(sizex, sizey))
    nA = len(A)
    nDB = len(DB)
    avg = nDB*A.astype(dtype=np.float16)
    for i in range(nDB):
        for j in range(DB.shape[1]):
            t = tuple(DB[i, j])
            avg[t] += 1.
            avg[t[::-1]] += 1.
    if nDB > 0:
        avg /= nDB

    # set colormap and colorbar
    cmap = plt.get_cmap('RdBu')
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=1, vmax=2))
    sm.set_array([])
    divider = make_axes_locatable(axs)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    ticks = (1.00, 1.25, 1.5, 1.75, 2.00);
    plt.colorbar(sm, label='averaged bond order', cax=cax, ticks=ticks)

    if annotate:
        for i in range(nA):
            axs.annotate(i, (V[i, 0], V[i, 1]))
    for i in range(nA):
        idx = np.transpose(np.nonzero(A[i]))
        for j in range(len(idx)):
            color = cmap(float(avg[i, idx[j]]-1.))
            lwidth = 9.*avg[i, idx[j]] - 8. # remormalize to [1,10]
            axs.plot((V[i, 0], V[idx[j], 0]),
                     (V[i, 1], V[idx[j], 1]),
                     c=color, ls='-', lw=lwidth)
    axs.set_xlim(min(V[:, 0])-2., max(V[:, 0])+2.)
    axs.set_ylim(min(V[:, 1])-2., max(V[:, 1])+2.)
    axs.set_xlabel('x [Ang]')
    axs.set_ylabel('y [Ang]')
    axs.set_aspect('equal')
    axs.patch.set_facecolor("black")
    if figname:
        plt.savefig(figname, dpi=dpi)
        plt.clf()
        plt.close(fig)
    else:
        plt.show()


def viewTBBondOrder(V, TB, figname=None, sizex=5, sizey=5, dpi=150, annotate=False):
    """
    Visualize the bond order estimated from tight-binding approach.
    """
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    fig, axs = plt.subplots(figsize=(sizex, sizey))
    nTB = len(TB)

    # set colormap and colorbar
    cmap = plt.get_cmap('RdBu')
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=1, vmax=2))
    sm.set_array([])
    divider = make_axes_locatable(axs)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    ticks = (1.00, 1.25, 1.5, 1.75, 2.00);
    plt.colorbar(sm, label='tight-binding bond order', cax=cax, ticks=ticks)

    # renormalize TB values to [1,2] interval
    vmin = np.min(TB[np.nonzero(TB)])
    vmax = np.max(TB)
    coef = 1./(vmax - vmin)
    TB *= coef
    coef = (vmax - 2.*vmin)*coef
    TB[np.nonzero(TB)] += coef

    if annotate:
        for i in range(nTB):
            axs.annotate(i, (V[i, 0], V[i, 1]))
    for i in range(nTB):
        idx = np.transpose(np.nonzero(TB[i]))
        for j in range(len(idx)):
            color = cmap(float(TB[i, idx[j]]-1))
            lrenorm = 9.*TB[i, idx[j]] - 8. # remormalize to [1,10]
            axs.plot((V[i, 0], V[idx[j], 0]),
                     (V[i, 1], V[idx[j], 1]),
                     c=color, ls='-', lw=lrenorm)
    axs.set_xlim(min(V[:, 0])-2., max(V[:, 0])+2.)
    axs.set_ylim(min(V[:, 1])-2., max(V[:, 1])+2.)
    axs.set_xlabel('x [Ang]')
    axs.set_ylabel('y [Ang]')
    axs.set_aspect('equal')
    axs.patch.set_facecolor("black")
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
