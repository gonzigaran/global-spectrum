#!/usr/bin/env python
# -*- coding: utf8 -*-

from utils import extend_const_sys, extend_non_sol_sys, antichain


def all_global_kernels(A, sigma, all_solutions=True):
    """
    Dada un algebra y un conjunto sigma de congruencias de A, devuelve todos
    los nucleos globales relativos a sigma
    """
    n = len(sigma)
    print(n)
    H_old = [([], A.maxcon(), [])]
    H_new = H_old.copy()
    solutions = []
    for i in range(n):
        for (deltas, intersection, vectors) in H_old:
            if antichain(sigma, i, deltas):
                intersection_new = intersection & sigma[i]
                deltas_new = deltas.copy()
                deltas_new.append(i)
                vectors_new = extend_const_sys(sigma, intersection, deltas, i)
                for vector in vectors:
                    vectors_new = vectors_new + extend_non_sol_sys(sigma,
                                                                   deltas,
                                                                   i,
                                                                   vector)
                new_tuple = (deltas_new, intersection_new, vectors_new)
                if intersection_new == A.mincon() and vectors_new == []:
                    # Si vectors_new == [] entonces la descomposición es global
                    # para el cociente de intersection_new
                    if not all_solutions:
                        return deltas_new
                    solutions.append(deltas_new)
                H_new.append(new_tuple)
        H_old = H_new
        print(i, len(H_old))
    return solutions


if __name__ == "__main__":
    import doctest

    doctest.testmod()
