#!/usr/bin/env python
# -*- coding: utf8 -*-

from utils import extend_const_sys, extend_non_sol_sys


def delta(ConQ_A, alphas, atom):
    """

    """
    joins = [ConQ_A.join(atom, alpha) for alpha in alphas]
    result = joins[0]
    for join in joins:
        result = ConQ_A.meet(result, join)
    return result


def is_global_indecomposable_atomics(ConQ_A, return_atomics=False):
    """
    Dada un algebra, decide si es globalmente indescomponible, usando el
    teorema de atomics

    >>> from folpy.examples.lattices import *
    >>> is_global_indecomposable_atomics(gen_chain(2).congruence_lattice())
    True
    >>> is_global_indecomposable_atomics(gen_chain(3).congruence_lattice())
    True
    >>> is_global_indecomposable_atomics(gen_chain(5).congruence_lattice())
    False
    >>> is_global_indecomposable_atomics(rhombus.congruence_lattice())
    False
    >>> is_global_indecomposable_atomics(M3.congruence_lattice())
    True
    >>> is_global_indecomposable_atomics(N5.congruence_lattice())
    True

    """
    A = ConQ_A.algebra
    thetas = ConQ_A.universe.copy()
    atoms = ConQ_A.atoms.copy()
    thetas.remove(A.mincon())
    n = len(thetas)
    pre2 = [[i, j]
            for i in range(n)
            for j in range(i + 1, n)
            if ConQ_A.meet(thetas[i], thetas[j]) == A.mincon()
            ]
    slots_old = [([i], []) for i in range(n)]  # Conjunto de la recursión,
    # contine pares de tuplas preatomicas junto a sus sistemas sin solución
    slots_new = slots_old.copy()
    for i in range(1, n):
        alpha = thetas[i]
        for (pre_atomic, vectors) in slots_old:
            # 0. Chequeo que las preatomicas queden ordenadas por el indice
            if i <= pre_atomic[-1]:
                continue

            # 1. Chequeo que con la nueva congruencia, la tupla sea preatomica
            is_new_pre_atomic = all([j, i] in pre2 for j in pre_atomic)
            if not is_new_pre_atomic:
                continue
            new_pre_atomic = pre_atomic + [i]

            # 2. Test de futuro
            has_future = True
            not_related_atoms = []
            atoms_copy = atoms.copy()
            while atoms_copy and has_future:
                atom = atoms_copy.pop()
                if ConQ_A.le(atom, alpha):
                    delta_i_atom = delta(ConQ_A,
                                         [thetas[i] for i in pre_atomic],
                                         atom)
                    has_future = ConQ_A.le(alpha, delta_i_atom)
                elif any(ConQ_A.le(atom, thetas[j]) for j in pre_atomic):
                    alpha_j_list = [j for j in pre_atomic
                                    if ConQ_A.le(atom, thetas[j])]
                    assert len(alpha_j_list) == 1
                    alpha_j = thetas[alpha_j_list[0]]
                    has_future = ConQ_A.le(alpha_j, ConQ_A.join(atom, alpha))
                else:
                    not_related_atoms.append(atom)

            # 3. Extender los vectores
            vectors_new = extend_const_sys(ConQ_A,
                                           thetas,
                                           A.mincon(),
                                           pre_atomic,
                                           i)
            for vector in vectors:
                vectors_new = vectors_new + \
                                extend_non_sol_sys(ConQ_A,
                                                   thetas,
                                                   pre_atomic,
                                                   i,
                                                   vector)

            # 4. Chequear A descomponible con test de presente
            if not vectors_new:  # Todo sistema para new_pre_atomic tiene
                # solucion
                is_decomposable = True
                # Test de presente
                while not_related_atoms and is_decomposable:
                    atom = not_related_atoms.pop()
                    delta_i_atom = delta(ConQ_A,
                                         [thetas[i] for i in new_pre_atomic],
                                         atom)
                    if all(not ConQ_A.le(thetas[j], delta_i_atom)
                            for j in new_pre_atomic):  # Si new_pre_atomic
                        # no es atomica
                        is_decomposable = False
                if is_decomposable:  # new_pre_atomic es atomica y todo sistema
                    # tiene solucion
                    if return_atomics:
                        [thetas[i] for i in new_pre_atomic]
                    else:
                        return False

            # 5. Agregar lote
            slots_new.append((new_pre_atomic, vectors_new))

        slots_old = slots_new.copy()

    return True
