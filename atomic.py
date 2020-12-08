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


def is_global_indecomposable_atomics(
                ConQ_A,
                return_atomics=False,
                verbose=False):
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
    mincon = A.mincon()
    thetas = ConQ_A.universe.copy()
    atoms = ConQ_A.atoms.copy()
    thetas.remove(A.mincon())
    thetas_enum = list(enumerate(thetas))
    pre2 = [[(i, theta), (j, delta)]
            for (i, theta) in thetas_enum
            for (j, delta) in thetas_enum[i:]
            if ConQ_A.meet(theta, delta) == mincon
            ]  # Conjunto con pares de congruencias disjuntas
    slots_old = [([(i, theta)], []) for (i, theta) in thetas_enum]  # Conjunto
    # de la recursión, contine pares de tuplas preatomicas junto a sus sistemas
    # sin solución
    slots_new = slots_old.copy()
    for (i, alpha) in thetas_enum:
        for (pre_atomic, vectors) in slots_old:
            # 0. Chequeo que las preatomicas queden ordenadas por el indice
            if i <= pre_atomic[-1][0]:
                continue

            # 1. Chequeo que con la nueva congruencia, la tupla sea preatomica
            is_new_pre_atomic = all([(j, delta), (i, alpha)] in pre2
                                    for (j, delta) in pre_atomic)
            if not is_new_pre_atomic:
                continue
            new_pre_atomic = pre_atomic + [(i, alpha)]

            # 2. Test de futuro
            has_future = True
            not_related_atoms = []
            atoms_copy = atoms.copy()
            while atoms_copy and has_future:
                atom = atoms_copy.pop()
                if ConQ_A.le(atom, alpha):
                    delta_i_atom = delta(ConQ_A,
                                         [theta for (j, theta) in pre_atomic],
                                         atom)
                    has_future = ConQ_A.le(alpha, delta_i_atom)
                elif any(ConQ_A.le(atom, theta) for (j, theta) in pre_atomic):
                    alpha_j_list = [(j, theta) for (j, theta) in pre_atomic
                                    if ConQ_A.le(atom, theta)]
                    assert len(alpha_j_list) == 1
                    alpha_j = alpha_j_list[0][1]
                    has_future = ConQ_A.le(alpha_j, ConQ_A.join(alpha, atom))
                else:
                    not_related_atoms.append(atom)

            # 3. Extender los vectores
            vectors_new = extend_const_sys(ConQ_A,
                                           [theta for (j, theta) in pre_atomic],
                                           alpha)
            for vector in vectors:
                vectors_new = vectors_new + \
                                extend_non_sol_sys(
                                    ConQ_A,
                                    [theta for (j, theta) in pre_atomic],
                                    alpha,
                                    vector
                                    )
            # 4. Chequear A descomponible con test de presente
            if not vectors_new:  # Todo sistema para new_pre_atomic tiene
                # solucion
                is_decomposable = True
                # Test de presente
                while not_related_atoms and is_decomposable:
                    atom = not_related_atoms.pop()
                    delta_i_atom = delta(
                        ConQ_A,
                        [theta for (j, theta) in new_pre_atomic],
                        atom
                        )
                    if all(not ConQ_A.le(theta, delta_i_atom)
                            for (j, theta) in new_pre_atomic):  # Si
                        # new_pre_atomic no es atomica
                        is_decomposable = False
                if is_decomposable:  # new_pre_atomic es atomica y todo sistema
                    # tiene solucion
                    if return_atomics:
                        return [theta for (j, theta) in new_pre_atomic]
                    else:
                        return False

            # 5. Agregar lote
            slots_new.append((new_pre_atomic, vectors_new))

        slots_old = slots_new.copy()

    return True
