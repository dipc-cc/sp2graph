"""

:mod:`sp2graph.visual`
===================================

Visualisation and ploting functions.

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


def viewKekule(V, A, DB, C=None, rad=None, figname=None,
               sizex=5, sizey=5, dpi=150, annotate=False):
    """
    Visualize a single Kekule representation with vertices
    coordinates 'V', adjacency matrix 'A' and double-bonds 'DB'.
    If provided, constrained bonds 'C' are shown with the usual Kekule
    representation and radicals are marked with a dot next to the vertex.
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

    # constrained double bonds
    if C:
        C = np.array(C, dtype=np.uint8)
        allC = C.flatten()
        # single bonds around all constrained vertices
        for i in range(len(allC)):
            ic = allC[i] # constrained vertex
            idx = np.transpose(np.nonzero(A[ic]))
            for j in range(len(idx)):
                axs.plot((V[ic, 0], V[idx[j], 0]),
                         (V[ic, 1], V[idx[j], 1]), c='#00FF37', ls='-', lw=1.5)
        # double bonds
        for i in range(len(C)):
            par = lau.parallel(V[C[i, 0]], V[C[i, 1]])
            axs.plot((par[0][0], par[1][0]),
                     (par[0][1], par[1][1]), c='y', ls='-', lw=1.5)

    # radicals
    if rad:
        for i in range(len(rad)):
            idx = np.transpose(np.nonzero(A[rad[i]]))
            radmk = lau.ptOrtho(V[idx[0]][0], V[rad[i]], V[idx[1]][0])
            axs.scatter(radmk[0], radmk[1], s=15, c='k', marker='o')
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


def viewKekuleGrid(V, A, DB, C=None, rad=None, figname=None,
                   sizex=5, sizey=5, dpi=150, annotate=False):
    """
    Visualize a single Kekule representation with vertices
    coordinates 'V', adjacency matrix 'A' and double-bonds 'DB'.
    If provided, constrained bonds 'C' are shown with the usual Kekule
    representation and radicals are marked with a dot next to the vertex.
    """
    nA = len(A)
    nKek = len(DB)
    DB = np.array(DB, dtype=np.uint8)
    nDB = DB.shape[1]

    # define the grid for ploting the figures
    molsize = max(V[:, 0]) - min(V[:, 0]) + 4.
    cols = int(40//molsize)
    rows = nKek/cols + nKek%cols

    fig = plt.figure()
    fig.set_size_inches(sizex, sizey)

    # set the radical markers (same for all Kekules)
    if rad:
        nrad = len(rad)
        radmk = np.zeros(shape=[nrad, nrad], dtype=np.float16)
        for i in range(len(rad)):
            idx = np.transpose(np.nonzero(A[rad[i]]))
            radmk[i] = lau.ptOrtho(V[idx[0]][0], V[rad[i]], V[idx[1]][0])

    # set constrained double bonds (same for all Kekules)
    if C:
        C = np.array(C, dtype=np.uint8)
        allC = C.flatten()
    else:
        C = np.empty(shape=[0, 0], dtype=np.uint8)
        allC = np.empty(shape=[0], dtype=np.uint8)

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
        # radicals
        if rad:
            plt.scatter(radmk[:, 0], radmk[:, 1], s=15, c='k', marker='o')

        # single bonds around all constrained vertices
        for i in range(len(allC)):
            ic = allC[i] # constrained vertex
            idx = np.transpose(np.nonzero(A[ic]))
            for j in range(len(idx)):
                plt.plot((V[ic, 0], V[idx[j], 0]),
                         (V[ic, 1], V[idx[j], 1]), c='#00FF37', ls='-', lw=1.5)
        # double bonds
        for i in range(len(C)):
            par = lau.parallel(V[C[i, 0]], V[C[i, 1]])
            plt.plot((par[0][0], par[1][0]),
                     (par[0][1], par[1][1]), c='y', ls='-', lw=1.5)

        if annotate:
            for i in range(nA):
                plt.annotate(i, (V[i, 0], V[i, 1]))
        """
        if idb%cols == 0:
            plt.ylabel('y [Ang]')
        else:
            plt.ylabel('')
            plt.tick_params(labelleft='off')
        plt.xlim(min(V[:, 0])-2., max(V[:, 0])+2.)
        plt.ylim(min(V[:, 1])-2., max(V[:, 1])+2.)
        plt.xlabel('x [Ang]')
        """
        plt.box(False)
        plt.axis('off')
        plt.title('#%i (%i)'%(idb+1, nKek))
        plt.gca().set_aspect(1)

    if figname:
        plt.savefig(figname, dpi=dpi)
        plt.clf()
        plt.close(fig)
    else:
        plt.show()


def viewBondOrderAverage(V, A, DB, C=None, rad=None, figname=None,
                         sizex=5, sizey=5, dpi=150, annotate=False):
    """
    Visualize a single Kekule representation with vertices
    coordinates 'V', adjacency matrix 'A' and double-bonds 'DB'.
    If provided, constrained bonds 'C' are shown with the usual Kekule
    representation and radicals are marked with a dot next to the vertex.
    """
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    fig, axs = plt.subplots(figsize=(sizex, sizey))
    nA = len(A)
    nDB = len(DB)
    avg = nDB*A.astype(dtype=np.float16)
    DB = np.array(DB, dtype=np.uint8)

    if len(DB.flatten()) == 0:
        print('WARNING: no Kekule strucutre to average!')
        return
    for i in range(nDB):
        for j in range(DB.shape[1]):
            t = tuple(DB[i, j])
            avg[t] += 1.
            avg[t[::-1]] += 1.
    if nDB > 0:
        avg /= nDB

    # constrained double bonds
    if C:
        C = np.array(C, dtype=np.uint8)
        allC = C.flatten()
    else:
        C = np.empty(shape=[0, 0], dtype=np.uint8)
        allC = np.empty(shape=[0], dtype=np.uint8)

    # set colormap and colorbar
    #cmap = plt.get_cmap('RdBu')
    cmap = plt.get_cmap('autumn')
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=1, vmax=2))
    sm.set_array([])
    divider = make_axes_locatable(axs)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    ticks = (1.00, 1.25, 1.5, 1.75, 2.00);
    plt.colorbar(sm, label='Pauling bond order', cax=cax, ticks=ticks)

    # single bonds around all constrained vertices
    for i in range(len(allC)):
        ic = allC[i] # constrained vertex
        idx = np.transpose(np.nonzero(A[ic]))
        for j in range(len(idx)):
            axs.plot((V[ic, 0], V[idx[j], 0]),
                     (V[ic, 1], V[idx[j], 1]), c='g', ls='-', lw=1.5)
    # constrained double bonds
    for i in range(len(C)):
        par = lau.parallel(V[C[i, 0]], V[C[i, 1]])
        axs.plot((par[0][0], par[1][0]),
                 (par[0][1], par[1][1]), c='y', ls='-', lw=1.5)
    # plot averaged of unconstrained only
    Auc = np.delete(A, allC, axis=0)
    for i in range(len(allC)):
        Auc[:, i] = 0
    # get indexes of unconstrained vertices
    allV = np.arange(len(A))
    unc = np.where(np.isin(allV, allC, invert=True))[0]
    for i in range(len(Auc)):
        idx = np.transpose(np.nonzero(Auc[i]))
        iunc = unc[i] # unconstrained vertex
        for j in range(len(idx)):
            color = cmap(float(avg[iunc, idx[j]]-1.))
            lwidth = 9.*avg[iunc, idx[j]] - 8. # remormalize to [1,10]
            axs.plot((V[iunc, 0], V[idx[j], 0]),
                     (V[iunc, 1], V[idx[j], 1]),
                     c=color, ls='-', lw=lwidth)
    # radicals
    if rad:
        for i in range(len(rad)):
            idx = np.transpose(np.nonzero(A[rad[i]]))
            radmk = lau.ptOrtho(V[idx[0]][0], V[rad[i]], V[idx[1]][0])
            axs.scatter(radmk[0], radmk[1], s=30, c='y', marker='o')
    if annotate:
        # Write Pauling bond orders
        for i in range(nA):
            for j in np.where(A[i, :]==1)[0]:
                axs.annotate('%.2f'%avg[i, j], ((V[i, 0]+V[j, 0])/2, (V[i, 1]+V[j, 1])/2), color='w')

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


def viewTBBondOrder(V, BO, figname=None,
                    sizex=5, sizey=5, dpi=150, annotate=False):
    """
    Visualize the bond order estimated from tight-binding approach.
    """
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    fig, axs = plt.subplots(figsize=(sizex, sizey))
    nBO = len(BO)

    # set colormap and colorbar
    #cmap = plt.get_cmap('RdBu')
    cmap = plt.get_cmap('autumn')
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=1, vmax=2))
    sm.set_array([])
    divider = make_axes_locatable(axs)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    ticks = (1.00, 1.25, 1.5, 1.75, 2.00);
    plt.colorbar(sm, label='Huckel bond order', cax=cax, ticks=ticks)

    for i in range(nBO):
        idx = np.transpose(np.nonzero(BO[i]))
        for j in range(len(idx)):
            color = cmap(float(BO[i, idx[j]]-1))
            lrenorm = 9.*BO[i, idx[j]] - 8. # remormalize to [1,10]
            axs.plot((V[i, 0], V[idx[j], 0]),
                     (V[i, 1], V[idx[j], 1]),
                     c=color, ls='-', lw=lrenorm)
    if annotate:
        # Write Pauling bond orders
        for i in range(len(BO)):
            for j in np.where(BO[i, :] > 0)[0]:
                if j > i:
                    axs.annotate('%.2f'%BO[i, j], ((V[i, 0]+V[j, 0])/2, (V[i, 1]+V[j, 1])/2), color='w')

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
