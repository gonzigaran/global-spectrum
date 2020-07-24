#!/usr/bin/env python
# -*- coding: utf8 -*-

from folpy.semantics.congruences import sup_proj


def minimals(sigma):
    """
    Dado un conjunto de congruencias `sigma` devuelve las minimales
    """
    result = [tita for tita in sigma if all(
        [not(tita > delta) for delta in sigma]
        )]
    return result


def antichain(sigma, i, deltas):
    """
    Dado un conjuntos de congruencias `delta` y una congruencia nueva `i`,
    decide si `i` es incomparable con las congruencias de `delta`
    """
    for j in deltas:
        if sigma[i] <= sigma[j] or sigma[j] <= sigma[i]:
            return False
    return True


def gen_roots(elements, delta):
    """
    Dado un conjunto de elementos `elements` y una congruencia `delta`,
    devuelve los representantes en `delta` para cada elemento en `elements`
    """
    result = []
    for x in elements:
        root = delta.root(x)
        if root not in result:
            result.append(root)
    return result


def extend_const_sys(sigma, intersection, deltas, i):
    """
    Dado un conjunto de congruencias y una congruencia nueva, genera sistemas
    sin solucion con los primeros elementos constantes
    """
    output = []
    n = len(deltas)
    if n == 0:
        return output
    gammas = []
    for k in range(n):
        gammas.append(sup_proj(sigma,
                               sigma[i],
                               sigma[deltas[k]]))
    algebra = sigma[i].algebra
    for a in algebra.universe:
        blocks = [gammas[k].block(a) for k in range(n)]
        elements = frozenset.intersection(*blocks)
        elements_roots = gen_roots(elements, sigma[i])
        for z in elements_roots:
            blocks = [sigma[deltas[k]].block(a) for k in range(n)]
            blocks += [sigma[i].block(z)]
            solutions = frozenset.intersection(*blocks)
            if solutions == frozenset({}):
                a_tuple = n * [a]
                output.append(a_tuple + [z])
    return output


def extend_non_sol_sys(sigma, deltas, i, vector):
    """
    Dado un conjunto de sismtemas de congruencias sin solucion y una
    congruencia nueva, extiende el sistema a nuevos sistemas sin solucion
    """
    output = []
    n = len(deltas)
    blocks = [(sup_proj(sigma,
                        sigma[i],
                        sigma[deltas[k]])).block(vector[k]) for k in range(n)]
    elements = frozenset.intersection(*blocks)
    elements_roots = gen_roots(elements, sigma[i])
    for z in elements_roots:
        output.append(vector + [z])
    return output


if __name__ == "__main__":
    import doctest

    doctest.testmod()
