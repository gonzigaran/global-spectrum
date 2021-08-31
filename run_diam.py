#!/usr/bin/env python
# -*- coding: utf8 -*-

from datetime import datetime
from functools import reduce
import logging
from itertools import combinations, product

from folpy.utils.parser.parser import Parser

from globalspectrum import is_global_indecomposable
# from globalkernel import all_global_kernels
from atomic import is_global_indecomposable_atomics
# from minion_globalspectrum import is_global_spectrum_minion_relation_n



def maximal_substructures(
        model,
        supermodel=None,
        filter_isos=True,
        filter_subdirect=False):
    if not supermodel:
        supermodel = model
    if filter_subdirect:
        assert len(model.factors) > 1
    universe = model.universe.copy()
    result = []
    result_compl = []
    for i in range(1, len(model)):
        has_subs_of_this_len = False
        for subset in combinations(universe, i):
            subset = set(subset)
            if any(x.issubset(subset) for x in result_compl):
                continue
            has_subs_of_this_len = True
            possible_subuniverse = [x for x in universe if x not in subset]
            if is_subuniverse(possible_subuniverse, model):
                substructure = supermodel.restrict(possible_subuniverse)
                result_compl.append(subset)
                are_isos = any(substructure.is_isomorphic(x) for x in result)
                if filter_isos and are_isos:
                    continue
                result.append(substructure)
                yield substructure
                yield from maximal_substructures(
                                substructure,
                                supermodel=supermodel,
                                filter_isos=True,
                                filter_subdirect=False)
        if not has_subs_of_this_len:
            break
    return result

def is_subuniverse(possible_subuniverse, model):
    for x in product(possible_subuniverse, repeat=2):
        if model.join(*x) not in possible_subuniverse:
            return False
        if model.meet(*x) not in possible_subuniverse:
            return False
    return True






def check_isos(sub, subs):
    for s in subs:
        logging.debug(s.universe)
        if sub.is_isomorphic_graph(s):
            logging.debug("iso")
            return True
    logging.debug("no iso")
    return False


def gen_subdirect_sublattices(lattices, verbose=False):
    start_time = datetime.now()
    L = reduce((lambda x, y: x * y), lattices)
    subs = []
    j = 0
    i = 0
    for sub in maximal_substructures(L):
        j = j + 1
        if j % 1 == 0:
            logging.debug('Substructure Nº: %s (Time: %s)', j, (datetime.now() - start_time))
        is_subdirect = sub.is_subdirect()
        logging.debug("subdirect: %s", is_subdirect)
        is_iso = sub.is_isomorphic_graph(L)
        logging.debug("iso: %s", is_iso)
        if is_subdirect and not is_iso:
            logging.debug("cumple")
            if not check_isos(sub, subs):
                subs.append(sub)
                i = i + 1
                if i % 1 == 0:
                    logging.debug('Subdirect Product Nº: %s', i)
    if verbose:
        logging.info("Cantidad de sublattices subdirectos: %s", len(subs))
        logging.info("Tiempo sublattices: %s", (datetime.now() - start_time))

    return subs


if __name__ == "__main__":
    logging.basicConfig(
        filename='out_run_diam2.log',
        format='%(asctime)s - %(levelname)s: %(message)s',
        level=logging.DEBUG
        )

    start_time = datetime.now()
    C2 = Parser("examples/lattices/2chain.model").parse()
    M3 = Parser("examples/lattices/M3.model").parse()
    DiamDiam = M3 * M3
    Diam2 = M3 * C2

    logging.info("Tiempo carga de modelos: %s", (datetime.now() - start_time))
    logging.info("Tiempo total: %s", (datetime.now() - start_time))

    logging.info("--------------------------")
    logging.info("Diam*2")
    logging.info("--------------------------")
    subs_Diam2 = gen_subdirect_sublattices([M3, C2], verbose=True)

    i = 0
    for lat in subs_Diam2:
        lat.to_file('output/diam_variety/Diam2_%s.model' % i)
        # lat.draw()
        i += 1

    logging.info("--------------------------")
    logging.info("Diam*Diam")
    logging.info("--------------------------")
    subs_DiamDiam = gen_subdirect_sublattices([M3, M3], verbose=True)

    i = 0
    for lat in subs_DiamDiam:
        lat.to_file('output/diam_variety/DiamDiam_%s.model' % i)
        i += 1

    i = 0
    for lat in subs_Diam2:
        con_lat = lat.congruence_lattice()
        logging.info("Is %s globally indecomposable? %s" %
                     (i, is_global_indecomposable_atomics(con_lat)))
        logging.info("Is %s globally indecomposable? %s" %
                     (i, is_global_indecomposable(lat,
                                                  congruence_lattice=con_lat)))
        i += 1

    logging.info("FIN")

    # continuos = subs_Diam2[4].continous()
    # con = continuos[0].principal_congruence(1, 5)
    # print(continuos)
    # print(con)
    # con_lat = subs_Diam2[4].congruence_lattice()
    # con_lat.draw()
    # is_global_indecomposable_atomics(con_lat)


