#!/usr/bin/env python
# -*- coding: utf8 -*-


def minimals(sigma):
    """
    Dado un conjunto de congruencias `sigma` devuelve las minimales
    """
    result = [tita for tita in sigma if all(
        [not(tita > delta) for delta in sigma]
        )]
    return result


def antichain(deltas, tita):
    """
    Dado un conjuntos de congruencias `deltas` y una congruencia nueva `tita`,
    decide si `tita` es incomparable con las congruencias de `deltas`
    """
    for delta in deltas:
        if tita <= delta or delta <= tita:
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


def extend_const_sys(projective, deltas, tita):
    """
    Dado un conjunto de congruencias y una congruencia nueva, genera sistemas
    sin solucion con los primeros elementos constantes
    """
    gammas = []
    for delta in deltas:
        gammas.append(projective.join(tita, delta))
    universe = tita.algebra.universe.copy()
    output = []
    n = len(deltas)
    if n == 0:
        return output
    while universe:
        a = universe.pop()
        blocks = [gamma.block(a) for gamma in gammas]
        elements = frozenset.intersection(*blocks)
        blocks = [delta.block(a) for delta in deltas]
        solutions = frozenset.intersection(*blocks)
        elements_roots = gen_roots(elements - solutions, tita)
        for z in elements_roots:
            if solutions.intersection(tita.block(z)) == frozenset({}):
                a_tuple = n * [a]
                output.append(a_tuple + [z])
        universe = [x for x in universe if x not in solutions]
    return output


def extend_non_sol_sys(projective, deltas, tita, vector):
    """
    Dado un conjunto de sismtemas de congruencias sin solucion y una
    congruencia nueva, extiende el sistema a nuevos sistemas sin solucion
    """
    output = []
    n = len(deltas)
    blocks = [(projective.join(
                    tita,
                    deltas[k]
                            )).block(vector[k]) for k in range(n)]
    elements = frozenset.intersection(*blocks)
    elements_roots = gen_roots(elements, tita)
    for z in elements_roots:
        output.append(vector + [z])
    return output


if __name__ == "__main__":
    import doctest

    doctest.testmod()
