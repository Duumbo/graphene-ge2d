#!/usr/bin/env python
# -*- coding-UFT-8 -*-
"""
Script __main__.py
Porte d'entrée pour le programme graphene_ge2d
"""
import sys
from struc_bande import write_matrix, compute_eigenvalues, generate_graph


def run_matrix_gen():
    return write_matrix.__main__()


def run_eigenvalues_gen():
    return compute_eigenvalues.__main__()


def run_graph_gen():
    return generate_graph.__main__()


def __main__():
    switch = sys.argv[1]
    if switch == "-m":
        return run_matrix_gen()
    elif switch == "-e":
        return run_eigenvalues_gen()
    elif switch == "-g":
        return run_graph_gen()
    else:
        raise ValueError("Invalid command: " + switch)


if __name__ == "__main__":
    __main__()
