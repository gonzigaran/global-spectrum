#!/usr/bin/env python
# -*- coding: utf8 -*-

from functools import lru_cache

from utils import extend_const_sys, extend_non_sol_sys


#@lru_cache(maxsize=64)
def delta(ConQ_A, thetas, atoms, ix_tuple, i_atom):
    """

    """
    atom = atoms[i_atom]
    alphas = [thetas[i] for i in ix_tuple]
    joins = [ConQ_A.join(atom, alpha) for alpha in alphas]
    result = joins[0]
    for join in joins:
        result = ConQ_A.meet(result, join)
    return result


def is_global_indecomposable_atomics(ConQ_A):
    """
    Dada un algebra, decide si es globalmente indescomponible, usando el
    teorema de atomics

    >>> from folpy.examples.lattices import *
    >>> is_global_indecomposable_atomics(gen_chain(5).congruence_lattice())
    True

    """
    A = ConQ_A.algebra
    thetas = ConQ_A.universe.copy()
    atoms = ConQ_A.atoms.copy()
    thetas.remove(A.mincon())
    n = len(thetas)
    n_atoms = len(atoms)
    pre2 = [[i, j] 
                for i in range(n) 
                for j in range(i,n) 
                if (ConQ_A.meet(thetas[i], thetas[j]) == A.mincon()
                and thetas[i] != thetas[j])
            ]
    slots_old = [([i],[]) for i in range(n)]
    slots_new = slots_old.copy()
    for i in range(1, n):
        alpha = thetas[i]
        for (pre_atomic, vectors) in slots_old:
            # 0. Chequeo que las preatomicas queden ordenadas por el indice
            if i <= pre_atomic[-1]:
                continue

            # 1. Chequeo que con la nueva congruencia, la tupla sea preatomica 
            is_new_pre_atomic = all([j,i] in pre2 for j in pre_atomic)
            if not is_new_pre_atomic:
                continue
            new_pre_atomic = pre_atomic + [i]

            # 2. Test de futuro
            has_future = True
            not_related_atoms = []
            for i_atom in range(n_atoms):
                if ConQ_A.lt(atoms[i_atom], alpha):
                    delta_i_atom = delta(ConQ_A, thetas, atoms, pre_atomic, i_atom)
                    has_future = has_future and ConQ_A.lt(alpha, delta_i_atom)
                elif any(ConQ_A.lt(atoms[i_atom], thetas[x]) for x in pre_atomic):
                    alpha_j_list = [x for x in pre_atomic if ConQ_A.lt(atoms[i_atom], thetas[x])]
                    assert len(alpha_j_list) == 1
                    alpha_j = thetas[alpha_j_list[0]]
                    has_future = has_future and ConQ_A.lt(alpha_j, ConQ_A.join(atoms[i_atom], alpha))
                else:
                    not_related_atoms.append(i_atom)
            if not has_future:
                continue

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
            if not vectors_new:
                is_decomposable = True
                # Test de presente
                for i_atom in not_related_atoms:
                    delta_i_atom = delta(ConQ_A,
                                         thetas,
                                         atoms,
                                         new_pre_atomic,
                                         i_atom)
                    if all(not ConQ_A.lt(thetas[x], delta_i_atom) for x in new_pre_atomic):
                        is_decomposable = False
                        break
                if is_decomposable:
                    return False

            # 5. Agregar lote
            slots_new.append((new_pre_atomic, vectors_new))
        
        slots_old = slots_new.copy()
    
    return True