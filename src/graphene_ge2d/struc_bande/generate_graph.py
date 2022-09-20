#!/usr/bin/env python
# -*- coding-UFT-8 -*-
"""
Diag de la matrice pour la structure de bande
"""
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import sparse
from numba import jit
from .matrice_coef import gen_gs, get_matrice
from ._constantes import (PATH, NDIAG, SIZE_MATRICES, GAMMA, POINT_M, POINT_K)


def joli_graphique(fig, ax, *args):
    """
    Mise en forme du graphique
    fig: l'objet figure
    ax: l'axe actif
    """
    chemin_norm = args[0]
    position_m = chemin_norm[PATH[0] - 1]
    position_k = chemin_norm[PATH[0] + PATH[1] - 1]
    ax.set_xticks(
            [0., position_m, position_k, chemin_norm[-1]],
            [r"$\Gamma$", r"$M$", r"$K$", r"$\Gamma$"]
    )
    ax.set_title(f"Structure de bande pour {NDIAG=}, {SIZE_MATRICES=}")
    return fig, ax


def __main__():
    # Initialisation du graphique
    fig, ax = plt.subplots()

    # Load data
    list_eigvals = np.load("raw_data/list_eigvals.npy")
    chemin_norm = np.load("raw_data/chemin_norm.npy")

    for i in range(8):
        """
        Zone de production du graphique.
        Position temporaire dans le code, à déplacer dans un autre
        script.
        """
        # print(list_eigvals)
        ax.plot(chemin_norm, list_eigvals[:, i], "o", markersize=2)

    # Beauté du graphique
    fig, ax = joli_graphique(fig, ax, *[chemin_norm])

    plt.show()


if __name__ == "__main__":
    __main__()
