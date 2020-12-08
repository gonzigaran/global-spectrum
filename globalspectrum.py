#!/usr/bin/env python
# -*- coding: utf8 -*-

from folpy.semantics.lattices import Projective

from utils import extend_const_sys, extend_non_sol_sys, minimals
from globalkernel import all_global_kernels


def is_global_spectrum(A, sigma):
    """
    Dada un algebra `A` y un conjunto `sigma` de congruencias de `A`, decide si
    `sigma` es un espectro global de `A`

    >>> from folpy.examples.lattices import *
    >>> sigma = rhombus.congruences().copy()
    >>> sigma.remove(rhombus.mincon())
    >>> sigma.remove(rhombus.maxcon())
    >>> is_global_spectrum(rhombus, sigma)
    True
    """
    sigma_m = minimals(sigma)
    deltas = []
    vectors = []
    intersection = A.maxcon()
    projective = Projective(sigma)
    for tita in sigma_m:
        vectors_new = extend_const_sys(projective,
                                       deltas,
                                       tita)
        for vector in vectors:
            vectors_new = vectors_new + extend_non_sol_sys(projective,
                                                           deltas,
                                                           tita,
                                                           vector)
        deltas.append(tita)
        intersection = intersection & tita
        vectors = vectors_new
    if vectors == []:
        return True
    else:
        return False


def is_global_indecomposable(A, congruence_lattice=None, verbose=False):
    """
    Dada un algebra, decide si es globalmente indescomponible

    >>> from folpy.examples.lattices import *
    >>> is_global_indecomposable(gen_chain(2))
    True
    >>> is_global_indecomposable(gen_chain(3))
    True
    >>> is_global_indecomposable(gen_chain(4))
    False
    >>> is_global_indecomposable(gen_chain(5))
    False
    >>> is_global_indecomposable(rhombus)
    False
    """
    if congruence_lattice:
        sigma = congruence_lattice.universe.copy()
    else:
        sigma = A.congruences().copy()
    sigma.remove(A.mincon())
    return all_global_kernels(A,
                              sigma,
                              all_solutions=False,
                              sigma_full=True,
                              verbose=verbose,
                              projective=congruence_lattice) == []


if __name__ == "__main__":
    import doctest

    doctest.testmod()
