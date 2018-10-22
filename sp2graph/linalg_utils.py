"""

:mod:`sp2graph.linalg_utils`
===================================

Simple linear algebra functions.

.. currentmodule:: sp2graph.linalg_utils

"""

import numpy as np

__all__ = ['unitA', 'closeV', 'parallel', 'ptOrtho', 'findClosest']


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


def parallel(pt1, pt2):
    """
    Returns a pair of points in R^2 belonging to a line
    parallel to the one defined by the provided points.
    """
    theta = angleA(pt2-pt1, [1., 0.])
    parpt1 = (pt1[0] + 0.2*np.cos(theta), pt1[1] + 0.2*np.sin(theta))
    parpt2 = (pt2[0] - 0.2*np.cos(theta), pt2[1] - 0.2*np.sin(theta))
    theta = theta - np.pi/2
    parpt1 = (parpt1[0] - 0.25*np.cos(theta), parpt1[1] - 0.25*np.sin(theta))
    parpt2 = (parpt2[0] - 0.25*np.cos(theta), parpt2[1] - 0.25*np.sin(theta))
    return parpt1, parpt2


def ptOrtho(v1, v2, v3):
    """
    Return a point slightly distant from 'v2' and belonging
    to the line orthogonal to (\vec{v2-v1}+\vec{v3-v2})
    and that contains 'v2'.
    """
    v12 = v2 - v1
    v23 = v3 - v2
    theta = angleA(v12, v23)
    ortho = (v12 - v23) * np.sin(0.5*theta)
    ortho = unitA(ortho)
    return (v2 + 0.35*ortho)


def findClosest(array, value):
    """
    Returns the element index from 'array'
    that is closer to 'value'.
    """
    array = np.asarray(array)
    idx = np.searchsorted(array, value)
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
