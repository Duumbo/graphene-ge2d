#!/usr/bin/env python
# -*- coding-UFT-8 -*-
"""
Diag de la matrice pour la structure de bande
"""
from tqdm import tqdm
import numpy as np
from scipy import sparse
from scipy import linalg
from . import _constantes as cons  # NDIAG


def __main__():
    chemin = np.load("raw_data/chemin.npy")
    matrices = np.load("raw_data/diagonals.npy")

    list_eigvals = np.empty((chemin.shape[0], 7))

    for i, k in tqdm(enumerate(chemin)):
        # Calcul des valeurs propres avec la bibliothèque ARPACK
        eigval = linalg.eigh(
                matrices[:, :, i],
                subset_by_index=[0, 6],
                eigvals_only=True
        )
        # Tableau de sortie des valeurs propres
        list_eigvals[i, :] = eigval
    np.save("raw_data/list_eigvals.npy", list_eigvals)


if __name__ == "__main__":
    __main__()
