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
        logging.debug(s.universe)
        if sub.is_isomorphic(s):
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
    for sub in L.substructures(
            filter_isos=True,
            filter_subdirect=False,
            proper=True):
        j = j + 1
        if j % 1 == 0:
            logging.debug(
                'Substructure Nº: %s (Time: %s)',
                j,
                (datetime.now() - start_time))
        is_subdirect = sub.is_subdirect()
        logging.debug("subdirect: %s", is_subdirect)
        is_iso = sub.is_isomorphic(L)
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
        # lat.draw()
        i += 1

    logging.info("--------------------------")
    logging.info("Globalmente Indescomponibles")
    logging.info("--------------------------")

    i = 0
    logging.info("ver globalmente indescomponible para subestructuras de M3xC2")
    for lat in subs_Diam2:
        con_lat = lat.congruence_lattice()
        logging.info("Is %s globally indecomposable? %s" %
                     (i, is_global_indecomposable_atomics(con_lat)))
        logging.info("Is %s globally indecomposable? %s" %
                     (i, is_global_indecomposable(lat,
                                                  congruence_lattice=con_lat)))
        i += 1

    i = 0
    logging.info("ver globalmente indescomponible para subestructuras de M3xM3")
    for lat in subs_DiamDiam:
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
