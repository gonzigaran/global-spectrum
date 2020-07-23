#!/usr/bin/env python
# -*- coding: utf8 -*-

from utils import extend_const_sys, extend_non_sol_sys, minimals


def is_global_spectrum(A, sigma):
    """
    Dada un algebra `A` y un conjunto `sigma` de congruencias de `A`, decide si
    `sigma` es un espectro global de `A`
    """
    sigma_m = minimals(sigma)
    deltas = []
    vectors = []
    intersection = A.maxcon()
    for i in range(len(sigma_m)):
        vectors_new = extend_const_sys(sigma, intersection, deltas, i)
        for vector in vectors:
            vectors_new = vectors_new + extend_non_sol_sys(sigma,
                                                           deltas,
                                                           i,
                                                           vector)
        deltas.append(i)
        intersection = intersection & sigma[i]
        vectors = vectors_new
    if vectors == []:
        return True
    else:
        return False


if __name__ == "__main__":
    import doctest

    doctest.testmod()
