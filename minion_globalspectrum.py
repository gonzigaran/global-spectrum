#!/usr/bin/env python
# -*- coding: utf8 -*-

from itertools import combinations, product

from folpy.semantics.lattices import Projective
from folpy.utils.minion import MinionSol

from utils import minimals


def gen_non_sol_tuples(sigma_m, projective):
    universe = projective.algebra.universe
    for xs in product(universe, repeat=len(sigma_m)):
        if all(sigma_m[i].is_root(x) for i, x in enumerate(xs)):
            blocks = [sigma_m[i].block(x) for i, x in enumerate(xs)]
            intersection = frozenset.intersection(*blocks)
            if intersection == frozenset({}):
                yield xs


def is_global_spectrum_minion_relation_n(A, sigma, block=None):
    """
    Dada un algebra `A` y un conjunto `sigma` de congruencias de `A`, decide si
    `sigma` es un espectro global de `A`

    >>> from folpy.examples.lattices import *
    >>> sigma = rhombus.congruences().copy()
    >>> sigma.remove(rhombus.mincon())
    >>> is_global_spectrum_minion_relation_n(rhombus, sigma)
    True
    """
    assert A.is_continous(), "El universo de A tiene que ser de la forma \
                              [0...n]"
    sigma_m = minimals(sigma)
    projective = Projective(sigma)
    relation = []
    for rel_tuple in gen_non_sol_tuples(sigma_m, projective):
        relation.append(rel_tuple)
    if len(relation) == 0:
        return True

    to_minion = "MINION 3\n\n"
    to_minion += "**VARIABLES**\n"
    to_minion += "DISCRETE x[%s]{0..%s}\n\n" % (len(sigma_m), len(A) - 1)
    to_minion += "**TUPLELIST**\n"
    for i, j in combinations(range(len(sigma_m)), r=2):
        to_minion += "J%sJ%s %s 2\n" % (i,
                                        j,
                                        len(projective.join(sigma_m[i],
                                                            sigma_m[j]
                                                            ).table()))
        for a, b in projective.join(sigma_m[i], sigma_m[j]).table():
            to_minion += "%s %s\n" % (a, b)
        to_minion += "\n"
    to_minion += "J%s %s %s\n" % (len(sigma_m) + 1,
                                  len(relation),
                                  len(sigma_m))
    for rel_tuple in relation:
        to_minion += " ".join([str(x) for x in rel_tuple]) + "\n"
    to_minion += "\n"
    to_minion += "**CONSTRAINTS**\n"
    for i, j in combinations(range(len(sigma_m)), r=2):
        to_minion += "table([x[%s],x[%s]],J%sJ%s)\n" % (i, j, i, j)
    for rel_tuple in relation:
        xs = ",".join(["x[" + str(x) + "]" for x in rel_tuple])
        to_minion += "table([%s],J%s)\n" % (xs, len(sigma_m) + 1)
    to_minion += "\n\n"
    to_minion += "**EOF**"
    return len(MinionSol(to_minion))


if __name__ == "__main__":
    import doctest

    doctest.testmod()
