#!/usr/bin/env python
# -*- coding: utf8 -*-

from folpy.semantics.classes import check_isos
from folpy.utils.parser.parser import Parser

from globalspectrum import is_global_indecomposable
from globalkernel import all_global_kernels
# from minion_globalspectrum import is_global_spectrum_minion_relation_n


if __name__ == "__main__":
    c2 = Parser("examples/semilattices/2chain.model").parse()
    c3 = Parser("examples/semilattices/3chain.model").parse()
    A = Parser("examples/semilattices/A.model").parse()
    helmet = Parser("examples/semilattices/helmet.model").parse()

    sigma = helmet.congruences_in([c2, c3, A])

    subalgebras = helmet.substructures()
    subs = []
    for emb, sub in subalgebras:
        if len(sub) != len(helmet):
            if not check_isos(sub, subs, sub.type):
                subs.append(sub.continous()[0])

    print("Is globally indecomposable? %s" % is_global_indecomposable(helmet))

    print("Total subalgebras: %s" % len(subs))

    print("Total subalgebras that are globally indecomposable: %s\
        " % sum([is_global_indecomposable(s) for s in subs]))

    kernels = all_global_kernels(
        helmet,
        sigma,
        all_solutions=False,
        verbose=True)

    print("Global kernels: %s" % kernels)
