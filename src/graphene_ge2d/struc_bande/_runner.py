#!/usr/bin/env python
# -*- coding-UFT-8 -*-
"""
Script runner.py
    Interactive runner for the generation of graphs.
"""
from . import _constantes as cons
from . import (
        write_matrix, compute_eigenvalues, generate_graph, potentiel, coquilles
)

cons.update_constantes()


def run_all():
    write_matrix.__main__()
    compute_eigenvalues.__main__()
    generate_graph.update_struct_bande()
    #potentiel.update_pot_map()


def run_pot_map():
    potentiel.update_pot_map()


def run_struct_bande():
    write_matrix.__main__()
    compute_eigenvalues.__main__()
    generate_graph.update_struct_bande()


def __main__():
    write_matrix.__main__()
    compute_eigenvalues.__main__()
    generate_graph.__main__()
    #potentiel.__main__()
