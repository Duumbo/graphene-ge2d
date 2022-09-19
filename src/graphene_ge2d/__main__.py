#!/usr/bin/env python
# -*- coding-UFT-8 -*-
"""
Diag de la matrice pour la structure de bande
"""
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.linalg as linalg
from scipy import sparse
import scipy.special as sp
from numba import jit, types
import numba_scipy
from .matrice_coef import coeffs, coeffs_sim, gen_gs, get_matrice
from ._constantes import *


def joli_graphique(fig, ax, *args):
    """
    Mise en forme du graphique
    fig: l'objet figure
    ax: l'axe actif
    """
    chemin_norm = args[0]
    position_m = chemin_norm[PATH[0] - 1]
    position_k = chemin_norm[PATH[0] + PATH[1] - 1]
    ax.set_xticks([0., position_m, position_k, chemin_norm[-1]], [r"$\Gamma$", r"$M$", r"$K$", r"$\Gamma$"])
    ax.set_title(f"Structure de bande pour {NDIAG=}, {SIZE_MATRICES=}")
    return fig, ax


def tracer_chemin():

    # Chemin 1, c'est de Gamma à M
    partie_chemin_1 = np.linspace(GAMMA, POINT_M, PATH[0], endpoint = False)

    # Chemin 2, c'est de M à K
    partie_chemin_2 = np.linspace(POINT_M, POINT_K, PATH[1], endpoint = False)

    # Chemin 3, c'est de K à Gamma
    partie_chemin_3 = np.linspace(POINT_K, GAMMA, PATH[2], endpoint = False)

    # On combine les chemins entre-eux
    chemin = np.append(partie_chemin_1, partie_chemin_2, axis = 0)
    chemin = np.append(chemin, partie_chemin_3, axis = 0)

    # Running difference du chemin
    running_diff = [0]
    for i, elem in enumerate(chemin[:0:-1]):
        diff = np.linalg.norm(elem - chemin[-(i+2)])
        running_diff.append(running_diff[-1] + diff)
    return chemin, running_diff


def __main__():
    fig, ax = plt.subplots()

    chemin, chemin_norm = tracer_chemin()

    # Tableau où l'on place les valeurs propres obtenues
    list_eigvals = np.empty((chemin.shape[0], 10))
    vecteurs, size_matrix = gen_gs(SIZE_MATRICES)
    try:
        for i, k in tqdm(enumerate(chemin)):
            diagonals = get_matrice(size_matrix, k, NDIAG, vecteurs)
            matrix_to_diag = sparse.diags([diagonals[:, 0], diagonals[:-1, 1], diagonals[:-1, 2]], [0, 1, -1])
            if False:
                """
                Debug, pour afficher la heat map des valeurs absolues des matrices.
                """
                np.savetxt(f"{i}.txt", matrix_to_diag)
                sns.heatmap(np.abs(matrix_to_diag), ax = ax)
                plt.show()
                fig.savefig("heatmap_matrice_F.png")

            eigval = sparse.linalg.eigsh(matrix_to_diag, k = 10, return_eigenvectors = False, which="SM")
            if np.any(eigval > 1e90):
                print(i)
                break
            list_eigvals[i, :] = eigval
        for i in range(10):
            #print(list_eigvals)
            ax.plot(chemin_norm, list_eigvals[:, i], "o", markersize=2)

            # Beauté du graphique
            fig, ax = joli_graphique(fig, ax, *[chemin_norm])

        plt.show()
        np.savetxt("raw_data/list_eigvals.txt.gz", list_eigvals)
    except Exception as e:
        print("Error, safety save")
        np.savetxt("raw_data/safety_save.txt.gz", list_eigvals)
        print(e)


if __name__ == "__main__":
    __main__()

