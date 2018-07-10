"""

:mod:`sp2graph.geometry`
===================================

Function for reading carbon sp2 geometries and project into xy plane.

.. currentmodule:: sp2graph.geometry

"""

import numpy as np
import sp2graph.linalg_utils as lau

__all__ = ['readgeom']


def readgeom(ifile):
    """
    Read a carbon based sp2 geometry and rotate it to the xy plane.
    """
    V = readCs(ifile)
    V = xyprojA(V)
    return V


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
        if row[0] == 'C':
            V = np.append(V, [[row[1], row[2], row[3]]], axis=0)
    f.close()
    V = V.astype(np.float32)
    return V


def xyprojA(array):
    """ Returns the an equiarrayalent 2D structure rotated to the xy plane """
    u1 = lau.unitA(array[1]-array[0])
    u2 = array[2]-array[0]
    u3 = lau.unitA(np.cross(u1, u2))
    u2 = lau.unitA(np.cross(u3, u1))
    proj = np.linalg.inv([u1, u2, u3])
    array = np.matmul(array, proj)
    if np.max(array[:, 2])-np.min(array[:, 2]) > 0.01:
        raise ValueError('Structure is not planar!')
    return np.delete(array, 2, axis=1)


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
