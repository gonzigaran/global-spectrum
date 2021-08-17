#!/usr/bin/env python
# -*- coding: utf8 -*-

from datetime import datetime
from functools import reduce
import logging

from folpy.utils.parser.parser import Parser

from globalspectrum import is_global_indecomposable
# from globalkernel import all_global_kernels
from atomic import is_global_indecomposable_atomics
# from minion_globalspectrum import is_global_spectrum_minion_relation_n


def check_isos(sub, subs):
    for s in subs:
        if sub.is_isomorphic_graph(s):
            return True
    return False


def gen_subdirect_sublattices(lattices, verbose=False):
    start_time = datetime.now()
    L = reduce((lambda x, y: x * y), lattices)
    subs = []
    j = 0
    i = 0
    for emb, sub in L.substructures():
        j = j + 1
        if j % 500 == 0:
            logging.debug('%s (Time: %s)', j, (datetime.now() - start_time))
        if sub.is_subdirect() and not sub.is_isomorphic(L):
            if not check_isos(sub, subs):
                subs.append(sub)
                i = i + 1
                if i % 10 == 0:
                    logging.debug(i)
    if verbose:
        logging.info("Cantidad de sublattices subdirectos: %s", len(subs))
        logging.info("Tiempo sublattices: %s", (datetime.now() - start_time))

    return subs


if __name__ == "__main__":
    logging.basicConfig(
        filename='run_diam.log',
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

    # continuos = subs_Diam2[4].continous()
    # con = continuos[0].principal_congruence(1, 5)
    # print(continuos)
    # print(con)
    # con_lat = subs_Diam2[4].congruence_lattice()
    # con_lat.draw()
    # is_global_indecomposable_atomics(con_lat)
