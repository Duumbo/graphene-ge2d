#!/usr/bin/env python
# -*- coding-UFT-8 -*-
"""
Diag de la matrice pour la structure de bande
"""
import numpy as np
from numba import jit
from .matrice_coef import get_matrice
from . import _constantes as cons  # (PATH, NDIAG, SIZE_MATRICES, GAMMA, POINT_M, POINT_K)


def tracer_chemin():

    # Chemin 1, c'est de Gamma à M
    partie_chemin_1 = np.linspace(cons.GAMMA, cons.POINT_M, cons.PATH[0], endpoint=False)

    # Chemin 2, c'est de M à K
    partie_chemin_2 = np.linspace(cons.POINT_M, cons.POINT_K, cons.PATH[1], endpoint=False)

    # Chemin 3, c'est de K à Gamma
    partie_chemin_3 = np.linspace(cons.POINT_K, cons.GAMMA, cons.PATH[2], endpoint=False)

    # On combine les chemins entre-eux
    chemin = np.append(partie_chemin_1, partie_chemin_2, axis=0)
    chemin = np.append(chemin, partie_chemin_3, axis=0)

    # Running difference du chemin
    # Pour pondérer l'axe des x sur le graphique
    running_diff = [0]
    for i, elem in enumerate(chemin[:0:-1]):
        diff = np.linalg.norm(elem - chemin[-(i+2)])
        running_diff.append(running_diff[-1] + diff)
    return chemin, running_diff


#@jit(nopython=True)
def produire_diagonales_no_scipy(path_k, vects):
    size = cons.n_mat
    output = np.empty((size, size, len(path_k)), dtype=np.complex128)
    for i, k in enumerate(path_k):
        output[:, :, i] = get_matrice(k, size, vects)
    return output


def __main__():
    # Production du chemin en k
    chemin, chemin_norm = tracer_chemin()

    # Tableau contenant les vecteurs G
    vecteurs = np.load("raw_data/vecteurs.npy")
    #vecteurs_x = np.multiply(vecteurs[:, 0], cons.b1[:, np.newaxis])
    #vecteurs_y = np.multiply(vecteurs[:, 1], cons.b2[:, np.newaxis])
    # Tableau contenant les diagonales des matrices sparses à diagonaliser
    matrices = produire_diagonales_no_scipy(chemin, vecteurs)

    np.save("raw_data/chemin.npy", chemin)
    np.save("raw_data/chemin_norm.npy", chemin_norm)
    np.save("raw_data/diagonals.npy", matrices)


if __name__ == "__main__":
    __main__()
