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
    zero_value = cons.v_0 * cons.r_0_bar / (cons.eps * cons.eta)

    chemin = np.load("raw_data/chemin.npy")
    diagonals = np.load("raw_data/diagonals.npy")

    # Positions des diagonales dans la matrice sparse
    positions_diag = [
            diag
            if diag % 2 == 0
            else -diag
            for diag in range(cons.NDIAG)
    ]

    list_eigvals = np.empty((chemin.shape[0], 8))

    for i, k in tqdm(enumerate(chemin)):
        # Dépack du tableau diagonales en matrices sparses
        matrix_to_diag = sparse.diags(
                [
                    diagonals[:, diag, i]
                    for diag in range(cons.NDIAG)
                ],
                positions_diag
        ).toarray()
        matrix_to_diag[np.isclose(matrix_to_diag, 0, atol=0)] = zero_value

        # Calcul des valeurs propres avec la bibliothèque ARPACK
        eigval = linalg.eigh(
                matrix_to_diag,
                subset_by_index=[0, 7],
                eigvals_only=True
        )
        # Tableau de sortie des valeurs propres
        list_eigvals[i, :] = eigval
    np.save("raw_data/list_eigvals.npy", list_eigvals)


if __name__ == "__main__":
    __main__()
