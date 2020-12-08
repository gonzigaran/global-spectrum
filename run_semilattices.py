#!/usr/bin/env python
# -*- coding: utf8 -*-

from datetime import datetime

from folpy.semantics.classes import check_isos
from folpy.utils.parser.parser import Parser

from globalspectrum import is_global_indecomposable
from globalkernel import all_global_kernels
from atomic import is_global_indecomposable_atomics
# from minion_globalspectrum import is_global_spectrum_minion_relation_n


if __name__ == "__main__":
    start_time = datetime.now()
    c2 = Parser("examples/semilattices/2chain.model").parse()
    c3 = Parser("examples/semilattices/3chain.model").parse()
    A = Parser("examples/semilattices/A.model").parse()
    helmet = Parser("examples/semilattices/helmet.model").parse()

    print("--------------------------")
    print("Tiempo carga de modelos: %s" % (datetime.now() - start_time))
    print("Tiempo total: %s" % (datetime.now() - start_time))
    print("--------------------------")
    print("\n")

    start_sigma = datetime.now()
    sigma = helmet.congruences_in([c2, c3, A])

    print("--------------------------")
    print("Tiempo para crear sigma: %s" % (datetime.now() - start_sigma))
    print("Tiempo total: %s" % (datetime.now() - start_time))
    print("--------------------------")
    print("\n")

    subalgebras = helmet.substructures()
    subs = []
    for emb, sub in subalgebras:
        if len(sub) != len(helmet):
            if not check_isos(sub, subs, sub.type):
                subs.append(sub.continous()[0])

    start_con_helmet = datetime.now()
    con_helm = helmet.congruence_lattice()
    print("--------------------------")
    print("Tiempo congruence lattice helmet: %s" %
          (datetime.now() - start_con_helmet))
    print("Tiempo total: %s" % (datetime.now() - start_time))
    print("--------------------------")
    print("\n")

    start_helmet_atomic = datetime.now()
    print("Is globally indecomposable? %s" %
          is_global_indecomposable_atomics(con_helm))

    print("--------------------------")
    print("Tiempo global helmet atomic: %s" %
          (datetime.now() - start_helmet_atomic))
    print("Tiempo total: %s" % (datetime.now() - start_time))
    print("--------------------------")
    print("\n")

    start_helmet_kernel = datetime.now()
    print("Is globally indecomposable? %s" %
          is_global_indecomposable(helmet, congruence_lattice=con_helm))

    print("--------------------------")
    print("Tiempo global helmet kernels: %s" %
          (datetime.now() - start_helmet_kernel))
    print("Tiempo total: %s" % (datetime.now() - start_time))
    print("--------------------------")
    print("\n")

    print("Total subalgebras: %s" % len(subs))
    print("\n")

    start_helmet_atomic = datetime.now()
    print("Total subalgebras that are globally indecomposable atomic: %s" %
          sum([is_global_indecomposable_atomics(s.congruence_lattice())
               for s in subs]))
    print("--------------------------")
    print("Tiempo globals subalgebras atomic: %s" %
          (datetime.now() - start_helmet_atomic))
    print("Tiempo total: %s" % (datetime.now() - start_time))
    print("--------------------------")
    print("\n")

    start_helmet_kernel = datetime.now()
    print("Total subalgebras that are globally indecomposable kernels: %s" %
          sum([is_global_indecomposable(s) for s in subs]))
    print("--------------------------")
    print("Tiempo globals subalgebras kernels: %s" %
          (datetime.now() - start_helmet_kernel))
    print("Tiempo total: %s" % (datetime.now() - start_time))
    print("--------------------------")

    print("\n\n\n")
    print("Tiempo total: %s" % (datetime.now() - start_time))
    print("\n\n\n")
    kernels = all_global_kernels(
        helmet,
        sigma,
        all_solutions=False,
        verbose=True)

    print("Global kernels: %s" % kernels)
