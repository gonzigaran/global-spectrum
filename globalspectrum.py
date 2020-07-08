#!/usr/bin/env python
# -*- coding: utf8 -*-


def antichain(sigma, i, deltas):
    """
    Dado un conjuntos de congruencias, decide si es una anticadena o no
    """
    for j in deltas:
        if sigma[i] <= sigma[j] or sigma[j] <= sigma[i]:
            return False
    return True


def does_extend(sigma, deltas, x_tuple, i, z):
    assert len(deltas) == len(x_tuple)
    n = len(deltas)
    for j in range(n):
        join = sigma[deltas[j]] | sigma[i]
        if x_tuple[j] not in join.block(z):
            return False
    return True


def has_solution(sigma, deltas, x_tuple, i, z):
    n = len(deltas)
    inters_class = sigma[i].block(z)
    for j in range(n):
        inters_class = inters_class & sigma[deltas[j]].block(x_tuple[j])
    return inters_class != frozenset({})


def extend_const_sys(sigma, intersection, deltas, i):
    """
    Dado un conjunto de congruencias y una congruencia nueva, genera sistemas
    sin solucion con los primeros elementos constantes
    """
    output = []
    n = len(deltas)
    for a in intersection.roots():
        clase_inter = intersection.block(a)
        for z in sigma[i].roots():
            clase_i = sigma[i].block(z)
            if clase_inter & clase_i == frozenset({}):
                a_tuple = n * [a]
                if does_extend(sigma, deltas, a_tuple, i, z):
                    output.append(a_tuple + [z])
    return output


def extend_non_sol_sys(sigma, deltas, i, vector):
    """
    Dado un conjunto de sismtemas de congruencias sin solucion y una
    congruencia nueva, extiende el sistema a nuevos sistemas sin solucion
    """
    output = []
    for z in sigma[i].roots():
        if does_extend(sigma, deltas, vector, i, z):
            if not has_solution(sigma, deltas, vector, i, z):
                output.append(vector + [z])
    return output


def all_subspectra(A, sigma, all_solutions=True):
    """
    Dada un algebra y un conjunto sigma de congruencias de A, devuelve todos
    los subconjuntos de sigma que son espectros globales de A
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
                        return new_tuple
                    solutions.append(new_tuple)
                H_new.append(new_tuple)
        H_old = H_new
        print(i, len(H_old))
    return solutions


if __name__ == "__main__":
    import doctest

    doctest.testmod()
