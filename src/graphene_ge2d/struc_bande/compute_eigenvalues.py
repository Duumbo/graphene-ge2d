#!/usr/bin/env python
# -*- coding-UFT-8 -*-
"""
Diag de la matrice pour la structure de bande
"""
from tqdm import tqdm
import numpy as np
from scipy import sparse
from ._constantes import NDIAG


def __main__():
    chemin = np.load("raw_data/chemin.npy")
    diagonals = np.load("raw_data/diagonals.npy")

    # Positions des diagonales dans la matrice sparse
    positions_diag = [
            diag
            if diag % 2 == 0
            else -diag
            for diag in range(NDIAG)
    ]

    list_eigvals = np.empty((chemin.shape[0], 8))

    for i, k in tqdm(enumerate(chemin)):
        # Dépack du tableau diagonales en matrices sparses
        matrix_to_diag = sparse.diags(
                [
                    diagonals[:, diag, i]
                    for diag in range(NDIAG)
                ],
                positions_diag
        )

        # Calcul des valeurs propres avec la bibliothèque ARPACK
        eigval = sparse.linalg.eigsh(
                matrix_to_diag,
                k=8,
                return_eigenvectors=False,
                which="SM"
        )
        # Tableau de sortie des valeurs propres
        list_eigvals[i, :] = eigval
    np.save("raw_data/list_eigvals.npy", list_eigvals)


if __name__ == "__main__":
    __main__()
