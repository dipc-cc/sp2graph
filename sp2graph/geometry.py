"""

:mod:`sp2graph.geometry`
===================================

Function for reading carbon sp2 geometries and project into xy plane.

.. currentmodule:: sp2graph.geometry

"""

import numpy as np
from numpy import linalg as LA
import sp2graph.linalg_utils as lau

__all__ = ['readgeom', 'readLattice', 'readCs']


def readgeom(ifile, rotate=True):
    """
    Read a carbon based sp2 geometry and rotate it to the xy plane.
    """
    V = readCs(ifile)
    L = readLattice(ifile)
    if rotate:
        return xyprojA(V, L)
    return V, L


def readLattice(ifile):
    r"""
    If present, reads the lattice vectors. For example, if the
    xyz file contains the following second line:

    `Lattice="79.4 0.0 0.0 0.0 20.3 0.0 0.0 0.0 25.0"`

    then the lattice vectors are assigned as:

    .. math::

        \begin{array}{crrrc}
        \vec{L}_1 = [& 79.4 &  0.0 &  0.0 & ]\\
        \vec{L}_2 = [&  0.0 & 20.3 &  0.0 & ]\\
        \vec{L}_3 = [&  0.0 &  0.0 & 25.0 & ]
        \end{array}

    """
    L = None
    f = open(ifile)
    line = f.readline()
    line = f.readline() # analyze second line
    f.close()
    line = line.lower()
    line = line.replace('"', ' ')
    splitline = line.split()
    for i, s in enumerate(splitline):
        if "lattice" in s or "cell" in s:
            L = splitline[i+1:i+10]
            L = np.array(L, dtype=np.float16)
            L = L.reshape((3, 3))
            break

    return L


def readCs(ifile):
    """
    Read the 'n' carbon atoms from a xyz file and return
    a (3,n) numpy array containing their coordinates.
    """
    f = open(ifile)
    next(f)
    next(f)
    V = np.empty(shape=[0, 3])
    for l in f:
        row = l.split()
        if len(row) > 0: # ignore empty lines
            if row[0] == 'C':
                V = np.append(V, [[row[1], row[2], row[3]]], axis=0)
    f.close()
    V = V.astype(np.float32)
    return V


def xyprojA(xyz, Lxyz):
    """
    Returns an equivalent 2D structure rotated to the `xy` plane.
    If the lattice vectors are provided, apply the same rotation.
    """
    if len(xyz) < 3:
        print ('WARNING: I don\'t know how to rotate it, sorry!\n')
        return xyz, Lxyz

    # define the rotation
    u1 = lau.unitA(xyz[1]-xyz[0])
    u2 = xyz[2]-xyz[0]
    u3 = lau.unitA(np.cross(u1, u2))
    u2 = lau.unitA(np.cross(u3, u1))
    rot = LA.inv([u1, u2, u3])

    # rotate to the `xy` plane and discard `z` component
    rxyz = np.matmul(xyz, rot)
    if np.max(rxyz[:, 2])-np.min(rxyz[:, 2]) > 0.01:
        raise ValueError('Structure is not planar!')

    # lattice vectors
    rLxyz = Lxyz
    if np.any(Lxyz):
        # apply the same rotation on the lattice parameters
        for i, Li in enumerate(Lxyz):
            rLxyz[i] = np.matmul(Li, rot)

    return rxyz, rLxyz


def tests():
    """ tests """
    def treadCs():
        from os import listdir
        from os.path import isfile, join
        gpath = '../examples/geometries/'
        geoms = [f for f in listdir(gpath) if isfile(join(gpath, f))]
        for g in geoms:
            V = readCs(gpath+g)
            if not V.size > 0:
                raise ValueError('Error at \'readCs\'!')
    treadCs()

if __name__ == '__main__':
    tests()
