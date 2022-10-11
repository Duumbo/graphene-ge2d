#!/usr/bin/env python
# -*- coding-UFT-8 -*-
"""
Script __main__.py
Porte d'entrée pour le programme graphene_ge2d
"""
import matplotlib.pyplot as plt
import sys
from struc_bande import (
        write_matrix, compute_eigenvalues, generate_graph, potentiel, coquilles, _runner, sliders
)


def run_interactive():
    sliders._activate()
    _runner.__main__()
    plt.show()


def run_coquille_gen():
    return coquilles.__main__()


def run_matrix_gen():
    return write_matrix.__main__()


def run_potentiel_map():
    return potentiel.__main__()


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
    elif switch == "-p":
        return run_potentiel_map()
    elif switch == "-c":
        return run_coquille_gen()
    elif switch == "-i":
        run_interactive()
    else:
        raise ValueError("Invalid command: " + switch)


if __name__ == "__main__":
    __main__()
