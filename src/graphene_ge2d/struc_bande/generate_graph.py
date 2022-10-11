#!/usr/bin/env python
# -*- coding-UFT-8 -*-
"""
Diag de la matrice pour la structure de bande
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from . import _constantes as cons  # (PATH, NDIAG, SIZE_MATRICES)


this = sys.modules[__name__]


def joli_graphique(fig, ax, *args):
    """
    Mise en forme du graphique
    fig: l'objet figure
    ax: l'axe actif
    """
    chemin_norm = args[0]
    position_m = chemin_norm[cons.PATH[0] - 1]
    position_k = chemin_norm[cons.PATH[0] + cons.PATH[1] - 1]
    ax.set_xticks(
            [0., position_m, position_k, chemin_norm[-1]],
            [r"$\Gamma$", r"$M$", r"$K$", r"$\Gamma$"]
    )
    ax.set_ylabel(r"Énergie (Unités de $e_0$)")
    ax.set_title(f"Structure de bande")
    ax.set_ylim(-0.05, 8)
    return fig, ax


def __main__():
    # Initialisation du graphique
    fig, ax = plt.subplots()

    # Load data
    list_eigvals = np.loadtxt("raw_data/list_eigvals.txt")
    chemin_norm = np.load("raw_data/chemin_norm.npy")

    # Lines to store the graphs
    this.lines = []

    for i in range(list_eigvals.shape[1]):
        """
        Zone de production du graphique.
        Position temporaire dans le code, à déplacer dans un autre
        script.
        """
        # print(list_eigvals)
        this.lines.append(ax.plot(chemin_norm, list_eigvals[:, i], "o", markersize=1)[0])

    # Beauté du graphique
    fig, ax = joli_graphique(fig, ax, *[chemin_norm])

    def update_struct_bande():
        # Reload data
        list_eigvals = np.loadtxt("raw_data/list_eigvals.txt")
        chemin_norm = np.load("raw_data/chemin_norm.npy")

        for i, l in enumerate(this.lines):
            l.set_data(chemin_norm, list_eigvals[:, i])

        fig.canvas.draw_idle()

    this.update_struct_bande = update_struct_bande

    #fig.savefig("Images/structure_de_bande.png")


if __name__ == "__main__":
    __main__()
