"""

:mod:`sp2graph.linalg_utils`
===================================

Simple linear algebra functions.

.. currentmodule:: sp2graph.linalg_utils

"""

import numpy as np

__all__ = ['unitA', 'closeV']


def unitA(array):
    """ Returns the unit array """
    return array / np.linalg.norm(array)


def angleA(a1, a2):
    """ Returns the angle in radians between arrays 'a1' and 'a2' """
    a1_u = unitA(a1)
    a2_u = unitA(a2)
    return np.arccos(np.clip(np.dot(a1_u, a2_u), -1.0, 1.0))


def closeV(iv, V, rad):
    """ Returns an array off 1st neighbours indexes """
    idx = np.empty(shape=[0], dtype=np.int32)
    for i in range(iv):
        if np.linalg.norm(V[i]-V[iv]) < rad:
            idx = np.append(idx, [[i]])
    for i in range(iv+1, len(V)):
        if np.linalg.norm(V[i]-V[iv]) < rad:
            idx = np.append(idx, [[i]])
    return idx


def tests():
    """ tests """
    def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
        return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

    def tunitA():
        i = 0
        while i < 10:
            a = np.random.rand(3)
            ua = unitA(a)
            if not isclose(np.linalg.norm(ua), 1.0):
                raise ValueError('Error at \'tunitA\'!')
            i += 1
    tunitA()

if __name__ == '__main__':
    tests()
