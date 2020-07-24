#!/usr/bin/env python
# -*- coding: utf8 -*-

from itertools import combinations

from folpy.semantics.congruences import sup_proj
from folpy.utils.minion import MinionSol

from utils import minimals


def is_global_spectrum_minion(A, sigma):
    """
    Dada un algebra `A` y un conjunto `sigma` de congruencias de `A`, decide si
    `sigma` es un espectro global de `A`

    >>> from folpy.examples.lattices import *
    >>> sigma = rhombus.congruences()
    >>> sigma.remove(rhombus.mincon())
    >>> is_global_spectrum_minion(rhombus, sigma)
    True
    """
    sigma_m = minimals(sigma)
    joins = dict()
    for i, j in combinations(range(len(sigma_m)), r=2):
        joins[(i, j)] = sup_proj(sigma, sigma_m[i], sigma_m[j])
    
    to_minion = "MINION 3\n\n"
    to_minion += "**VARIABLES**\n"
    to_minion += "DISCRETE x[%s]{0..%s}\n\n" % (len(sigma_m) + 1, len(A) - 1)
    to_minion += "**TUPLELIST**\n"
    for (i, j) in joins:
        to_minion += "J%sJ%s %s 2\n" % (i, j, len(joins[(i, j)].table()))
        for a, b in joins[(i, j)].table():
            to_minion += "%s %s\n" % (a, b)
        to_minion += "\n"
    for i in range(len(sigma_m)):
        to_minion += "J%sJ%s %s 2\n" % (i, i, len(sigma_m[i].table()))
        for a, b in sigma_m[i].table():
            to_minion += "%s %s\n" % (a, b)
        to_minion += "\n"
    to_minion += "**CONSTRAINTS**\n"
    for (i, j) in joins:
        to_minion += "table([x[%s],x[%s]],J%sJ%s)\n" % (i, j, i, j)
    for i in range(len(sigma_m)):
        to_minion += "table([x[%s],x[%s]],J%sJ%s)\n" % (i, len(sigma_m), i, i)
    to_minion += "\n\n"
    to_minion += "**EOF**"
    return len(MinionSol(to_minion)) == len(A) ** len(sigma_m)


if __name__ == "__main__":
    import doctest

    doctest.testmod()
