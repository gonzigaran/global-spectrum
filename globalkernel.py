#!/usr/bin/env python
# -*- coding: utf8 -*-

from datetime import datetime
from folpy.semantics.lattices import Projective

from utils import extend_const_sys, extend_non_sol_sys, antichain


def all_global_kernels(A,
                       sigma,
                       all_solutions=True,
                       sigma_full=False,
                       verbose=False,
                       projective=None):
    """
    Dada un algebra y un conjunto sigma de congruencias de A, devuelve todos
    los nucleos globales relativos a sigma
    """
    if verbose:
        start_time = datetime.now()
    n = len(sigma)
    H_old = [([], A.maxcon(), [])]
    H_new = H_old.copy()
    solutions = []
    if not projective:
        projective = Projective(sigma + [A.mincon()], full=sigma_full)
    if verbose:
        print("Time gen projective: %s" % (datetime.now() - start_time))
        print("--------------------------------------------------")
    for tita in sigma:
        if verbose:
            i = sigma.index(tita)
            print("delta nº: %s (%s total)" % (i, n))
            print(sigma[i])
            print("Iterative set volume: %s" % len(H_old))
            print("Time: %s" % (datetime.now() - start_time))
            print("--------------------------------------------------")
        for (deltas, intersection, vectors) in H_old:
            if antichain(deltas, tita):
                intersection_new = intersection & tita
                deltas_new = deltas.copy()
                deltas_new.append(tita)
                vectors_new = extend_const_sys(projective,
                                               deltas,
                                               tita)
                for vector in vectors:
                    vectors_new = vectors_new + \
                                  extend_non_sol_sys(projective,
                                                     deltas,
                                                     tita,
                                                     vector)
                new_tuple = (deltas_new, intersection_new, vectors_new)
                if intersection_new == A.mincon() and vectors_new == []:
                    # Si vectors_new == [] entonces la descomposición es global
                    # para el cociente de intersection_new
                    if not all_solutions:
                        if verbose:
                            total_time = (datetime.now() - start_time)
                            print("Total time: %s" % total_time)
                        return deltas_new
                    solutions.append(deltas_new)
                H_new.append(new_tuple)
        H_old = H_new.copy()
    if verbose:
        total_time = (datetime.now() - start_time)
        print("Total time: %s" % total_time)
    return solutions


if __name__ == "__main__":
    import doctest

    doctest.testmod()
